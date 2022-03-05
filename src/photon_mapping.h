#pragma once

#include "scene.h"
#include "pcg.h"
#include "photon_map.h"
#include "material.h"
#include "util.h"



class PhotonMapping {

private:

	// number of photons for making caustics photon map
	int num_caustic_photons;

	// number of photons used for radiance estimation by caustics photon map
	int num_caustic_estimation;

	// maximum depth of photon tracing, eye tracing
	int max_depth;

	PhotonMap caustic_photon_map;
	PhotonMap global_photon_map;

/*	inline Vector3 sampleCosineHemisphere(const Vector2& uv, float& pdf) {
		const float theta =
			0.5f * std::acos(std::clamp(1.0f - 2.0f * uv[0], -1.0f, 1.0f));
		const float phi = PI_MUL_2 * uv[1];
		const float cosTheta = std::cos(theta);
		pdf = PI_INV * cosTheta;
		return sphericalToCartesian(theta, phi);
	}

	Vector3 sampleDirection(const Light light, pcg32_state &rng,
						  float& pdf)  {
		Vector2 uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
		const Vector3 dir = sampleCosineHemisphere(uv, pdf);

		// transform direction from local to world
		return localToWorld(dir, surfInfo.dpdu, surfInfo.shadingNormal,
							surfInfo.dpdv);
	}*/

	bool is_diffuse(const Material& mat) {
		// return true if material is Lambertian or DisneyDiffuse
		return mat.index() == 0 || mat.index() == 3;
	}

	bool is_specular(const Material& mat) {
		// return true if material is RoughDielectric
		return mat.index() == 2;
	}

	PathVertex get_path_vertex(const Scene &scene, const PointAndNormal &point_on_light) {
		// TODO:  refactor this to generate path vertex without intersection. This will fail for complex geometry
		Vector3 org = make_zero_spectrum();
		Vector3 pertub_pt = point_on_light.position +  0.5 * normalize(org - point_on_light.position);
		Vector3 dir = normalize(pertub_pt - point_on_light.position);

		Ray ray{pertub_pt, dir, Real(0), infinity<Real>()};

		std::optional<PathVertex> pv_opt = intersect(scene, ray);
		assert(pv_opt);
		return *pv_opt;
	}

	// finds a light direction and initializes the throughput
	Ray sample_photon_ray(const Scene& scene,
						  Spectrum &throughput, pcg32_state &rng) {
		// sample a light
		Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
		Real light_w = next_pcg32_real<Real>(rng);
		Real shape_w = next_pcg32_real<Real>(rng);
		int light_id = sample_light(scene, light_w);
		const Light &light = scene.lights[light_id];

		// TODO:  enable light sources other than area light
		// fail if the light source is not DiffuseAreaLight
		assert(light.index() == 0);
		DiffuseAreaLight areaLight = std::get<DiffuseAreaLight>(light);
		const Shape &shape = scene.shapes[areaLight.shape_id];

		// sample point on light
		Vector3 ref_point = make_zero_spectrum();
		PointAndNormal point_on_light =
			sample_point_on_light(light, ref_point, light_uv, shape_w, scene);

		ShadingInfo shading_info = compute_shading_info(shape, get_path_vertex(scene, point_on_light));

		// sample light direction
		Vector3 light_dir = to_world(shading_info.shading_frame, sample_cos_hemisphere(light_uv));

		// TODO:  revisit light_dir_pdf
		const float theta =
			0.5f * std::acos(std::clamp(1.0f - 2.0f * (float) light_uv.x, -1.0f, 1.0f));
		Real light_dir_pdf = c_INVPI * std::cos(theta);

		Ray light_ray{point_on_light.position, light_dir, get_intersection_epsilon(scene), infinity<Real>()};

		Spectrum L = emission(light,
							  light_dir, // pointing outwards from light
							  0,
							  point_on_light, // dummy parameter for envmap
							  scene);
		Real light_pdf = light_pmf(scene, light_id) *
				  pdf_point_on_light(light, point_on_light, ref_point, scene) * light_dir_pdf;

		// TODO: revisit this
		throughput = (L / light_pdf) *
					 std::abs(dot(light_dir, shading_info.shading_frame.n));
//		throughput = L;

		return light_ray;
	}

	void build_caustic_photon_map(const Scene &scene) {
		std::vector<Photon> photons;

		for (int i = 0; i < num_caustic_photons; ++i) {
			pcg32_state rng = init_pcg32(i);
			Vector3 throughput = make_const_spectrum(1);
			Ray ray = sample_photon_ray(scene, throughput, rng);
//			print(ray.org, "light org");
//			print(ray.dir, "light dir");

			bool prev_specular = false;
			for (int k = 0; k < max_depth; ++k) {
				std::optional<PathVertex> _vertex = intersect(scene, ray);
				if (_vertex) {
					PathVertex vertex = *_vertex;
					// evaluate bsdf at this point
					const Material &mat = scene.materials[vertex.material_id];
					// break when hitting diffuse surface without previous specular
					if (!prev_specular && is_diffuse(mat)) {
						break;
					}
					if (prev_specular && is_diffuse(mat)) {
						photons.emplace_back(throughput, vertex.position, -ray.dir);
						break;
					}
					prev_specular = is_specular(mat);

					// russian roulette
					// TODO:  try using rr_depth instead of doing it for k > 0
					if (k > 0) {
						const Real rr_prob = fmin(max(throughput), 1.0);
						if (next_pcg32_real<Real>(rng) > rr_prob) {
							// Terminate the path
							break;
						}
						throughput /= rr_prob;
					}
					Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
					Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

					// TODO:  the TransportDirection is not being set as TO_VIEW
					std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, -ray.dir, vertex,
																			   scene.texture_pool, bsdf_rnd_param_uv,
																			   bsdf_rnd_param_w, TransportDirection::TO_VIEW);
					if (!bsdf_sample_) {
						// BSDF sampling failed. Abort the loop.
						break;
					}
					const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
					Vector3 dir = bsdf_sample.dir_out;
					Real dir_pdf = pdf_sample_bsdf(mat, -ray.dir, dir, vertex, scene.texture_pool, TransportDirection::TO_VIEW);
/*					Real G = fabs(dot(ray.dir, vertex.geometry_normal)) /
							 distance_squared(ray.org, vertex.position);
					// convert pdf from solid angle  to area measure
					// TODO:  this geometric term is causing issues?
					// dir_pdf *= G;*/
					Spectrum f = eval(mat, -ray.dir, dir, vertex, scene.texture_pool, TransportDirection::TO_VIEW);

					throughput *= f / dir_pdf;

					ray = Ray{vertex.position, dir, get_intersection_epsilon(scene), infinity<Real>()};

				} else {
					break;
				}
			}
		}

		printf("generated %lu caustic photons\n", photons.size());
		caustic_photon_map.set_photons(photons);
	}

	void build_global_photon_map(const Scene &scene) {

	}

public:
	PhotonMapping(const Scene &scene) {
		this->num_caustic_photons = scene.options.num_caustic_photons;
		this->num_caustic_estimation = scene.options.num_caustic_estimation;

		this->max_depth = scene.options.max_depth == -1 ? 10 : scene.options.max_depth;
	}

	void build_photon_map(const Scene &scene) {
		build_caustic_photon_map(scene);
//		for (int i = 0; i < caustic_photon_map.get_num_photons(); ++i) {
//			print(caustic_photon_map.get_ith_photon(i).power, std::to_string(i + 1) + " caustic photon power");
//			print(caustic_photon_map.get_ith_photon(i).position, std::to_string(i + 1) + " caustic photon position");
//		}
		build_global_photon_map(scene);
	}

	Vector3 compute_caustics_with_photon_map(const Scene& scene, const Vector3 &wo, const PathVertex &vertex) {
		// get nearby photons
		Real max_dist2 = 0;
		const std::vector<int> photon_indices =
			caustic_photon_map.queryKNearestPhotons(vertex.position,
												   num_caustic_estimation, max_dist2);

		Vector3 Lo;
		for (const int photon_idx : photon_indices) {
			const Photon& photon = caustic_photon_map.get_ith_photon(photon_idx);
			const Material &mat = scene.materials[vertex.material_id];
			const Vector3 f = eval(mat, wo, photon.dir_in, vertex, scene.texture_pool, TransportDirection::TO_LIGHT);
			Lo += f * photon.power;
		}
		if (!photon_indices.empty()) {
			Lo /= (num_caustic_photons * c_PI * max_dist2);
		}

		return Lo;
	}

	Vector3 est_radiance_recursively(const Ray& ray, const Scene& scene, int depth)  {
		if (depth >= max_depth) return make_zero_spectrum();

		int w = scene.camera.width, h = scene.camera.height;
		RayDifferential ray_diff = init_ray_differential(w, h);

		std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
		if (!vertex_) {
			return make_zero_spectrum();
		}
		PathVertex vertex = *vertex_;

		if (is_light(scene.shapes[vertex.shape_id])) {
			return emission(vertex, -ray.dir, scene);
		}

		const Material &mat = scene.materials[vertex.material_id];

		if (is_diffuse(mat)) {
			// currently only caustics
			// compute caustics illumination with caustics photon map
			const Vector3 Lc = compute_caustics_with_photon_map(scene, -ray.dir, vertex);

			return Lc;
		} else if (is_specular(mat)) {
			// implement this
			return make_zero_spectrum();
		}
		return make_zero_spectrum();

	}

	Spectrum estimate_radiance(const Scene &scene,
							int x, int y, /* pixel coordinates */
							pcg32_state &rng) {
		int w = scene.camera.width, h = scene.camera.height;
		Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
						   (y + next_pcg32_real<Real>(rng)) / h);
		Ray ray = sample_primary(scene.camera, screen_pos);

		return est_radiance_recursively(ray, scene, 0);

	}
};