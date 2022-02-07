#pragma once
#include "util.h"

#define INF infinity<Real>()
// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
	int w = scene.camera.width, h = scene.camera.height;
	Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
					   (y + next_pcg32_real<Real>(rng)) / h);
	Ray ray = sample_primary(scene.camera, screen_pos);
	auto ray_diff = RayDifferential{Real(0), Real(0)};

	std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

	if (vertex_) {
		PathVertex vertex = *vertex_;
		Medium medium = scene.media[vertex.exterior_medium_id];
		Real tfar = distance(ray.org, vertex.position);
		Spectrum transmittance = exp(-get_sigma_a(medium, vertex.position) * tfar);

		Spectrum Le = make_zero_spectrum();
		if (is_light(scene.shapes[vertex.shape_id]))  {
			Le =  emission(vertex, -ray.dir, scene);
		}
		return transmittance * Le;
	}
    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
	int w = scene.camera.width, h = scene.camera.height;
	Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
					   (y + next_pcg32_real<Real>(rng)) / h);
	Ray ray = sample_primary(scene.camera, screen_pos);
	auto ray_diff = RayDifferential{Real(0), Real(0)};

	std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
	Medium medium = scene.media[scene.camera.medium_id];
	// TODO: how would we get a generic point instead of hardcoding?
	Spectrum sigma_a = get_sigma_a(medium, make_zero_spectrum());
	Spectrum sigma_s = get_sigma_s(medium, make_zero_spectrum());
	Spectrum _sigma_t = sigma_a + sigma_s;
	Real sigma_t = _sigma_t.x;

	// calculate t_hit
	Real t_hit = infinity<Real>();

	if (vertex_) {
		PathVertex vertex = *vertex_;
		t_hit = distance(ray.org, vertex.position);
	}

	Real u = next_pcg32_real<Real>(rng);
	Real t = -log(1 - u)/sigma_t;

	if (t < t_hit) {
		Real trans_pdf = exp(-sigma_t * t) * sigma_t;
		Real transmittance = exp(-sigma_t * t);

		Vector3 p = ray.org + t * ray.dir;

		// equation 7
		// eval function
		Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
		Real light_w = next_pcg32_real<Real>(rng);
		Real shape_w = next_pcg32_real<Real>(rng);
		int light_id = sample_light(scene, light_w);
		const Light &light = scene.lights[light_id];
		PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
		Vector3 dir_in = -ray.dir;
		Vector3 dir_out = normalize(point_on_light.position - p);
		Spectrum rho = eval(get_phase_function(medium), dir_in, dir_out);

		Spectrum Le = emission(light, -dir_out, Real(0), point_on_light, scene);
		Vector3 p_prime = point_on_light.position;
		Vector3 n_p_prime = point_on_light.normal;

		Real dist = distance(p, p_prime);
		Real exp_val = exp(-sigma_t * dist);
		Ray vis_ray{p, dir_out, get_shadow_epsilon(scene), infinity<Real>()};

		std::optional<PathVertex> int_point_ = intersect(scene, vis_ray, ray_diff);
		Real visibility = 1;
		PathVertex int_point = *int_point_;
		if (fabs(distance(int_point.position, p) - dist) > 0.001) {
			visibility = 0;
		}

		Real G = max(-dot(dir_out, n_p_prime), Real(0)) / (dist * dist);

		Spectrum L_s1_estimate = rho * Le * exp_val * G * visibility;
		Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);

		assert(trans_pdf != 0);
		assert(L_s1_pdf != 0);

		return (transmittance/trans_pdf) * sigma_s * (L_s1_estimate/L_s1_pdf);
	} else {
		Real trans_pdf = exp(-sigma_t * t_hit);
		Real transmittance = exp(-sigma_t * t_hit);
		Spectrum Le = make_zero_spectrum();

		// we can assume that the ray is intersecting the object. since t_hit cannot be infinity
		PathVertex vertex = *vertex_;
		if (is_light(scene.shapes[vertex.shape_id]))  {
			Le =  emission(vertex, -ray.dir, scene);
		}
		assert(trans_pdf != 0);
		return (transmittance/trans_pdf) * Le;
	}
}

// TODO: what is the third argument here?
int update_medium_id(Ray ray, PathVertex vertex, int medium_id) {
	if (vertex.interior_medium_id != vertex.exterior_medium_id) {
		if (dot(ray.dir, vertex.geometry_normal) > 0) {
			medium_id = vertex.exterior_medium_id;
		} else {
			medium_id = vertex.interior_medium_id;
		}
	}
	return medium_id;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
	int w = scene.camera.width, h = scene.camera.height;
	Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
					   (y + next_pcg32_real<Real>(rng)) / h);
	Ray ray = sample_primary(scene.camera, screen_pos);
	auto ray_diff = RayDifferential{Real(0), Real(0)};

	int current_medium_id = scene.camera.medium_id;

	Spectrum current_path_throughput = make_const_spectrum(1);
	Spectrum radiance = make_zero_spectrum();
	int bounces = 0;
	int max_depth = scene.options.max_depth;


	while(true)
	{
		bool scatter = false;
/*		if (debug(x, y)) {
//			std::optional<PathVertex> val  = intersect(scene, ray, ray_diff);
			printf("New bounce started for bounce %d!!\n", bounces);
			debug(x, y);
		}*/
		std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
		Real t_hit = INF;
		PathVertex vertex;
		if (vertex_) {
			vertex = *vertex_;
			t_hit = distance(ray.org, vertex.position);
		}

		Medium medium = scene.media[current_medium_id];
		Real sigma_a = get_sigma_a(medium, make_zero_spectrum()).x;
		Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;
		Real sigma_t = sigma_a + sigma_s;
		Real trans_pdf = 1;
		Real transmittance = 1;
		if (current_medium_id != -1) {
			Real u = next_pcg32_real<Real>(rng);
			Real t = -log(1 - u)/sigma_t;
			if (t < t_hit) {
				trans_pdf = exp(-sigma_t * t) * sigma_t;
				transmittance = exp(-sigma_t * t);
				scatter = true;
				ray.org += t * ray.dir;
			} else {
				trans_pdf = exp(-sigma_t * t_hit);
				transmittance = exp(-sigma_t * t_hit);
				ray.org = vertex.position;
			}
		}
		current_path_throughput *= (transmittance/trans_pdf);
/*		if (debug(x, y)) {
//			std::optional<PathVertex> val  = intersect(scene, ray, ray_diff);
			print(current_path_throughput, "current through put during bounce " + std::to_string(bounces) +
			" using sigma_t " + std::to_string(sigma_t));
			debug(x, y);
		}*/
		bool hit_light = false;
		if (!scatter) {
			// set emission if there is an intersection with a light source
			if (t_hit != INF  && is_light(scene.shapes[vertex.shape_id])) {
				hit_light = true;
				radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
/*				if (debug(x, y)) {
					printf("Updating radiance in bounce %d\n", bounces);
				}*/
			}
		}

		if (bounces == max_depth - 1 and max_depth != -1) {
			break;
		}

		if (!scatter && t_hit != INF && vertex.material_id == -1) {
			// index-matching interface, skip through it
			current_medium_id = update_medium_id(ray, vertex, current_medium_id);
			bounces++;
			continue;
/*			if (debug(x, y)) {
				printf("Medium updated to %d in bounce %d\n", current_medium_id, bounces);
			}*/
		}

		if (scatter) {
			Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
			const Vector3 dir_in = -ray.dir;
/*			if (debug(x, y)) {
				printf("Medium used has value %d in bounce %d\n", current_medium_id, bounces);
			}*/
			std::optional<Vector3> next_dir_ = sample_phase_function(get_phase_function(scene.media[current_medium_id]),
																	dir_in, phase_rnd_param_uv);
			if (!next_dir_) {
				break;
			}

			const Vector3 &next_dir = *next_dir_;
			ray = Ray{ray.org, next_dir, get_intersection_epsilon(scene), infinity<Real>()};

			Spectrum phase_func_eval = eval(get_phase_function(scene.media[current_medium_id]), dir_in, next_dir);
			Real phase_func_pdf = pdf_sample_phase(get_phase_function(scene.media[current_medium_id]), dir_in, next_dir);
			current_path_throughput *= (phase_func_eval/phase_func_pdf) * sigma_s;

		} else if (hit_light) {
			break;
		}

		Real rr_prob = 1;
		if (bounces >= scene.options.rr_depth) {
			rr_prob = fmin(max(current_path_throughput), Real(0.95));
			if (next_pcg32_real<Real>(rng) > rr_prob) {
				// Terminate the path
				break;
			} else {
				current_path_throughput /= rr_prob;
			}
		}
		bounces += 1;
	}
    return radiance;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
