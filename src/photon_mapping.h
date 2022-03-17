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

	// number of photons for making caustics photon map
	int num_global_photons;

	// number of photons used for radiance estimation by caustics photon map
	int num_global_estimation;

	// maximum depth to estimate radiance by final gathering
	int final_gathering_depth;

	// maximum depth of photon tracing, eye tracing
	int max_depth;

	PhotonMap *caustic_photon_map;
	PhotonMap *global_photon_map;

	bool is_diffuse(const Material& mat) {
		// return true if material is Lambertian or DisneyDiffuse
		return mat.index() == 0 || mat.index() == 3;
	}

	bool is_specular(const Material& mat) {
		// return true if material is RoughDielectric
		return mat.index() == 2;
	}

	std::optional<PathVertex> get_path_vertex(const Scene &scene, const PointAndNormal &point_on_light) {
		//// TODO:  refactor this to generate path vertex without intersection. This might fail for complex geometry
		Vector3 org = make_zero_spectrum();
		Vector3 pertub_pt = point_on_light.position +  2 * get_intersection_epsilon(scene)  * point_on_light.normal;
		Vector3 dir = normalize(-point_on_light.normal);
		Ray ray{pertub_pt, dir, get_intersection_epsilon(scene), infinity<Real>()};
		std::optional<PathVertex> pv_opt = intersect(scene, ray);
		return pv_opt;
	}

	// finds a light direction and initializes the throughput
	std::optional<Ray> sample_photon_ray(const Scene& scene,
						  Spectrum &throughput, pcg32_state &rng) {
		// sample a light
		while (true)
		{
			Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
			Real light_w = next_pcg32_real<Real>(rng);
			Real shape_w = next_pcg32_real<Real>(rng);
			int light_id = sample_light(scene, light_w);
			const Light& light = scene.lights[light_id];

			// TODO:  enable light sources other than area light
			// fail if the light source is not DiffuseAreaLight
			assert(light.index() == 0);
			DiffuseAreaLight areaLight = std::get<DiffuseAreaLight>(light);
			const Shape& shape = scene.shapes[areaLight.shape_id];

			// sample point on light
			Vector3 ref_point = make_zero_spectrum();
			PointAndNormal point_on_light =
				sample_point_on_light(light, ref_point, light_uv, shape_w, scene);

			std::optional<PathVertex> pvrt = get_path_vertex(scene, point_on_light);
			if (!pvrt.has_value())
			{
				continue;
			}
			ShadingInfo shading_info = compute_shading_info(shape, *pvrt);

			// sample light direction
			Vector3 light_dir = to_world(shading_info.shading_frame, sample_cos_hemisphere(light_uv));

			// TODO:  revisit light_dir_pdf
			const float theta =
				0.5f * std::acos(std::clamp(1.0f - 2.0f * (float)light_uv.x, -1.0f, 1.0f));
			Real light_dir_pdf = c_INVPI * std::cos(theta);

			Ray light_ray{ point_on_light.position, light_dir, get_intersection_epsilon(scene), infinity<Real>() };

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
		
	}

	void build_caustic_photon_map(const Scene &scene) {
		std::vector<Photon> photons;

		for (int i = 0; i < num_caustic_photons; ++i) {
			pcg32_state rng = init_pcg32(i);
			Vector3 throughput = make_const_spectrum(1);
			std::optional<Ray> _ray = sample_photon_ray(scene, throughput, rng);
			//if (!_ray.has_value())
			//{
			//	continue;
			//}
			Ray ray = *_ray;
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
						photons.emplace_back(vertex.position, throughput, -ray.dir);
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
		caustic_photon_map = new PhotonMap(photons);
	}

	void build_global_photon_map(const Scene &scene) {
		std::vector<Photon> photons;

		for (int i = 0; i < num_global_photons; ++i) {
			pcg32_state rng = init_pcg32(i);
			Vector3 throughput = make_const_spectrum(1);
			std::optional<Ray> _ray = sample_photon_ray(scene, throughput, rng);
			//if (!_ray.has_value())
			//{
			//	continue;
			//}
			Ray ray = *_ray;
			//std::cout << i << std::endl;
			for (int k = 0; k < max_depth; ++k) {
				std::optional<PathVertex> _vertex = intersect(scene, ray);
				if (_vertex) {
					PathVertex vertex = *_vertex;
					// evaluate bsdf at this point
					const Material &mat = scene.materials[vertex.material_id];
					// break when hitting diffuse surface without previous specular
					if (is_diffuse(mat)) {
						photons.emplace_back(vertex.position, throughput, -ray.dir);
					}

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
					Spectrum f = eval(mat, -ray.dir, dir, vertex, scene.texture_pool, TransportDirection::TO_VIEW);

					throughput *= f / dir_pdf;

					ray = Ray{vertex.position, dir, get_intersection_epsilon(scene), infinity<Real>()};

				} else {
					break;
				}
			}
		}

		printf("generated %lu global photons\n", photons.size());
		global_photon_map = new PhotonMap(photons);
	}

public:
	PhotonMapping(const Scene &scene) {
		this->num_caustic_photons = scene.options.num_caustic_photons;
		this->num_caustic_estimation = scene.options.num_caustic_estimation;
		this->num_global_photons = scene.options.num_global_photons;
		this->num_global_estimation = scene.options.num_global_estimation;

		this->final_gathering_depth = scene.options.final_gathering_depth;

		this->max_depth = scene.options.max_depth;
	}

	void print_photon_map(PhotonMap photonMap) {
		for (int i = 0; i < photonMap.get_num_photons(); ++i) {
			print(photonMap.get_ith_photon(i).power, std::to_string(i + 1) + " photon power");
			print(photonMap.get_ith_photon(i).position, std::to_string(i + 1) + " photon position");
		}
	}
	void build_photon_map(const Scene &scene) {
		if (final_gathering_depth > 0) {
			build_caustic_photon_map(scene);
		}

		build_global_photon_map(scene);
	}

	Vector3 compute_indirect_illumination_recursive(const Scene& scene,
													const Vector3& dir_view,
													const PathVertex& vertex,
													pcg32_state &rng,
													int depth) {
		if (depth >= max_depth) return make_zero_spectrum();

		int w = scene.camera.width, h = scene.camera.height;
		Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
		Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

		Vector3 Li = make_zero_spectrum();
		// sample direction by BxDF
		const Material &mat = scene.materials[vertex.material_id];
		std::optional<BSDFSampleRecord> bsdf_sample_ =
			sample_bsdf(mat, dir_view, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
		if (!bsdf_sample_) {
			// BSDF sampling failed. Abort the loop.
			return make_zero_spectrum();
		}
		const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
		Vector3 dir_bsdf = bsdf_sample.dir_out;

		const Vector3 f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
		Real pdf_dir = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);

		// trace final gathering ray
		Ray ray_fg{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
		RayDifferential ray_diff = init_ray_differential(w, h);

		std::optional<PathVertex> next_vertex_ = intersect(scene, ray_fg, ray_diff);
		if (!next_vertex_) {
			return make_zero_spectrum();
		}
		PathVertex next_vertex = *next_vertex_;
		const Material &next_mat = scene.materials[next_vertex.material_id];
		if (is_diffuse(next_mat)) {
			Li += f *
				  compute_radiance_with_photon_map(scene, -ray_fg.dir, next_vertex) /
				  pdf_dir;
		}
		else if (is_specular(next_mat)) {
			Li += f *
				compute_indirect_illumination_recursive(
					scene, -ray_fg.dir, next_vertex, rng, depth + 1) /
				  pdf_dir;
		}

		return Li;
	}

	// compute indirect illumination with final gathering
	Vector3 compute_indirect_illumination(const Scene& scene, const Vector3& wo,
										  const PathVertex& vertex,
										  pcg32_state &rng, int depth) {
		// TODO: should the depth be zero?
		return compute_indirect_illumination_recursive(scene, wo, vertex, rng, depth);
	}

	// compute direct illumination with explicit light sampling(NEE)
	Vector3 compute_direct_illumination(const Scene& scene, const Vector3& dir_view,
										const PathVertex &vertex,
										pcg32_state &rng) {
		Vector3 Ld = make_zero_spectrum();

		// sample light
		Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
		Real light_w = next_pcg32_real<Real>(rng);
		Real shape_w = next_pcg32_real<Real>(rng);
		int light_id = sample_light(scene, light_w);
		const Light &light = scene.lights[light_id];

		// sample point on light
		PointAndNormal point_on_light =
			sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);

		Vector3 dir_light = normalize(point_on_light.position - vertex.position);

		Ray shadow_ray{vertex.position, dir_light,
					   get_shadow_epsilon(scene),
					   (1 - get_shadow_epsilon(scene)) *
					   distance(point_on_light.position, vertex.position)};
		if (!occluded(scene, shadow_ray)) {
			Real pdf_pos_light = light_pmf(scene, light_id) *
					  pdf_point_on_light(light, point_on_light, vertex.position, scene);
			// convert positional pdf to directional pdf
			const Vector3 wi = normalize(point_on_light.position - vertex.position);
			const Real r = length(point_on_light.position - vertex.position);
			const Real pdf_dir =
				pdf_pos_light * r * r / std::abs(dot(-wi, point_on_light.normal));

			Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);

			const Material &mat = scene.materials[vertex.material_id];
			Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);

			Ld = f * Le / pdf_dir;
		}

		return Ld;
	}

	// compute reflected radiance with global photon map
	Vector3 compute_radiance_with_photon_map(const Scene& scene, const Vector3& wo,
											 const PathVertex& vertex)  {
		// get nearby photons
		Real max_dist2;
		const std::vector<Photon> photons =
			global_photon_map->queryKNearestPhotons(vertex.position,
												 num_global_estimation, max_dist2);

		Vector3 Lo;
		for (Photon photon : photons) {
			const Material &mat = scene.materials[vertex.material_id];

			const Vector3 f = eval(mat, wo, photon.dir_in, vertex, scene.texture_pool, TransportDirection::TO_LIGHT);
			Lo += f * photon.power;
		}
		if (!photons.empty()) {
			Lo /= (num_global_photons * c_PI * max_dist2);
		}
		return Lo;
	}

	Vector3 compute_caustics_with_photon_map(const Scene& scene, const Vector3 &wo, const PathVertex &vertex) {
		// get nearby photons
		Real max_dist2 = 0;
		const std::vector<Photon> photons =
			caustic_photon_map->queryKNearestPhotons(vertex.position,
												   num_caustic_estimation, max_dist2);

		Vector3 Lo;
		for (Photon photon : photons) {
			const Material &mat = scene.materials[vertex.material_id];
			const Vector3 f = eval(mat, wo, photon.dir_in, vertex, scene.texture_pool, TransportDirection::TO_LIGHT);
			Lo += f * photon.power;
		}
		if (!photons.empty()) {
			Lo /= (num_caustic_photons * c_PI * max_dist2);
		}

		return Lo;
	}

	Vector3 est_radiance_recursively(const Ray& ray, const Scene& scene, int depth, pcg32_state &rng)  {
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

			if (depth >= final_gathering_depth) {
				return compute_radiance_with_photon_map(scene, -ray.dir, vertex);
			} else {
				// compute direct illumination by NEE
				const Vector3 Ld =
					compute_direct_illumination(scene, -ray.dir, vertex, rng);

				// compute caustics illumination with caustics photon map
				const Vector3 Lc = compute_caustics_with_photon_map(scene, -ray.dir, vertex);

				// compute indirect illumination with final gathering
				const Vector3 Li =
					compute_indirect_illumination(scene, -ray.dir, vertex, rng, depth);

				return Lc + Ld + Li;
			}

		} else if (is_specular(mat)) {
			Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
			Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);

			// sample direction by BxDF
			Vector3 dir_view = -ray.dir;
			const Material &mat = scene.materials[vertex.material_id];
			std::optional<BSDFSampleRecord> bsdf_sample_ =
				sample_bsdf(mat, dir_view, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
			if (!bsdf_sample_) {
				// BSDF sampling failed. Abort the loop.
				return make_zero_spectrum();
			}
			const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
			Vector3 dir_bsdf = bsdf_sample.dir_out;

			const Vector3 f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
			Real pdf_dir = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
			// recursively raytrace
			const Ray next_ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
			const Vector3 throughput = f / pdf_dir;

			return throughput *
				est_radiance_recursively(next_ray, scene, depth + 1, rng);
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
		Spectrum result = est_radiance_recursively(ray, scene, 0, rng);
		if (result.x != result.x) {
			return Vector3(0, 0, 0);
		} else {
			return result;
		}
	}
};