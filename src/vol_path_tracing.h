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


	while(true) {
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

		if (!scatter && t_hit != INF  && is_light(scene.shapes[vertex.shape_id])) {
			// set emission if there is an intersection with a light source
			radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
/*				if (debug(x, y)) {
					printf("Updating radiance in bounce %d\n", bounces);
				}*/
		}

		if (!scatter && t_hit != INF && vertex.material_id == -1) {
			// index-matching interface, skip through it
			ray = Ray{vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
			current_medium_id = update_medium_id(ray, vertex, current_medium_id);
			bounces++;
			continue;
/*			if (debug(x, y)) {
				printf("Medium updated to %d in bounce %d\n", current_medium_id, bounces);
			}*/
		}

		if (bounces == max_depth - 1 and max_depth != -1) {
			break;
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

		} else {
			break;
		}

		Real rr_prob = 1;
		if (bounces >= scene.options.rr_depth) {
			rr_prob = fmin(max(current_path_throughput), Real(0.95));
			if (next_pcg32_real<Real>(rng) > rr_prob) {
				// Terminate the path
				break;
			}
		}
		bounces += 1;
		current_path_throughput /= rr_prob;
	}
    return radiance;
}

bool intersects(std::optional<PathVertex> pv) {
	return pv.has_value();
}

Real get_sigma_t(const Scene &scene, int medium_id) {
	Medium medium = scene.media[medium_id];
	Real sigma_a = get_sigma_a(medium, make_zero_spectrum()).x;
	Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;
	return sigma_a + sigma_s;
}

Spectrum next_event_estimation(const Scene &scene, Vector3 p, int current_medium_id, PathVertex vertex, const Vector3 dir_in,
							   int bounces, pcg32_state &rng, bool surface, int x, int y) {
	// First, we sample a point on the light source.
	Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
	Real light_w = next_pcg32_real<Real>(rng);
	Real shape_w = next_pcg32_real<Real>(rng);
	int light_id = sample_light(scene, light_w);
	const Light &light = scene.lights[light_id];
	PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
	Vector3 p_prime = point_on_light.position;
	Vector3 dir_light = normalize(p_prime - p);
	Vector3 orig_point = p;

	auto ray_diff = RayDifferential{Real(0), Real(0)};

	// Compute transmittance to light. Skip through index-matching shapes.
	Spectrum T_light = make_const_spectrum(Real(1));
	int shadow_medium_id = current_medium_id;
	int shadow_bounces = 0;
	Real p_trans_dir = 1;

	int max_depth = scene.options.max_depth;
	while(true) {
		Ray shadow_ray{p, dir_light, get_shadow_epsilon(scene),
					   (1 - get_shadow_epsilon(scene)) *
					   distance(point_on_light.position, p)};
		std::optional<PathVertex> lvertex_ = intersect(scene, shadow_ray, ray_diff);
		Real next_t = distance(p, p_prime);

		PathVertex lvertex;
		if (intersects(lvertex_)) {
			lvertex = *lvertex_;
			next_t = distance(p, lvertex.position);
		}
		// we have taken the following if out of the intersects if block because for, the shadow ray need not pass
		// though an IMS all the time.
		// Account for the transmittance to next_t
		if(shadow_medium_id >= 0) {
			T_light *= exp(-get_sigma_t(scene, shadow_medium_id) * next_t);
			p_trans_dir *= exp(-get_sigma_t(scene, shadow_medium_id) * next_t);
		}

		if (!intersects(lvertex_)) {
			// Nothing is blocking, we’re done
			break;
		} else {
			// Something is blocking: is it an opaque surface?
			if (lvertex.material_id >= 0) {
				// we're blocked
				return make_zero_spectrum();
			}
			// otherwise, it’s an index-matching surface and we want to pass through -- this introduces one extra
			// connection vertex
			shadow_bounces++;
			if ((bounces + shadow_bounces) >= max_depth-1 && max_depth != -1) {
				// Reach the max no. of vertices
				return make_zero_spectrum();
			}
			shadow_medium_id = update_medium_id(shadow_ray, lvertex, shadow_medium_id);
			p += next_t * shadow_ray.dir;
		}
	}

	if (debug(x, y)) {
		debug(x, y);
	}

	if (max(T_light) > 0) {
		Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
					 distance_squared(point_on_light.position, orig_point);
		Real pdf_nee, pdf_dir;
		Spectrum contrib;
		if (surface) {
			const Material &mat = scene.materials[vertex.material_id];
			Spectrum f = eval(mat, dir_in, dir_light, vertex, scene.texture_pool);
			Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
			pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, orig_point, scene);
			contrib = T_light * G * f * L / pdf_nee;

			pdf_dir = pdf_sample_bsdf(mat, dir_in, dir_light, vertex, scene.texture_pool) * G * p_trans_dir ;
		} else {
			PhaseFunction phase_function = get_phase_function(scene.media[current_medium_id]);
			Spectrum f = eval(phase_function, dir_in, dir_light);
			Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
			pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, orig_point, scene);
			contrib = T_light * G * f * L / pdf_nee;

			pdf_dir = pdf_sample_phase(phase_function, dir_in, dir_light) * G * p_trans_dir;
		}

		// power heuristics
		Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_dir * pdf_dir);
		return w * contrib;
	}

	return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
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
	Real dir_pdf = 0; // in solid angle measure
	Spectrum nee_p_cache;
	Real multi_trans_pdf = 1;

	int max_depth = scene.options.max_depth;

	while(true) {
		bool scatter = false;
		std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
		Real t_hit = INF;
		PathVertex vertex;
		if (vertex_) {
			vertex = *vertex_;
			t_hit = distance(ray.org, vertex.position);
		}

		Real trans_pdf = 1;
		Real transmittance = 1;
		if (current_medium_id != -1) {
			Medium medium = scene.media[current_medium_id];
			Real sigma_a = get_sigma_a(medium, make_zero_spectrum()).x;
			Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;
			Real sigma_t = sigma_a + sigma_s;
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
			multi_trans_pdf *= trans_pdf;
		}
		current_path_throughput *= (transmittance/trans_pdf);


		// emission for phase function sampling
		if (!scatter && t_hit != INF  && is_light(scene.shapes[vertex.shape_id])) {

			// set emission if there is an intersection with a light source
			if (bounces == 0 ) {
				radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
			} else {
				// Need to account for next event estimation
				int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
				assert(light_id >= 0);
				const Light &light = scene.lights[light_id];
				PointAndNormal light_point{vertex.position, vertex.geometry_normal};

				Real pdf_nee = light_pmf(scene, light_id) *
							   pdf_point_on_light(light, light_point, nee_p_cache, scene);
				
				Real G = fmax(Real(0), -dot(ray.dir, vertex.geometry_normal)) /
					distance_squared(nee_p_cache, vertex.position);
				Real pdf_phase = dir_pdf * multi_trans_pdf * G;
				Real weight = (pdf_phase * pdf_phase) / (pdf_phase * pdf_phase + pdf_nee * pdf_nee);
				radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * weight;
			}
		}

		// index-matching interface, skip through it
		if (!scatter && t_hit != INF && vertex.material_id == -1) {
			current_medium_id = update_medium_id(ray, vertex, current_medium_id);
			ray = Ray{vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};

			bounces++;
			continue;
		}

		if (bounces >= max_depth - 1 and max_depth != -1) {
			break;
		}

		// radiance from nee will be added here
		if (scatter && current_medium_id != -1) {
			Medium medium = scene.media[current_medium_id];
			Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;

			const Vector3 dir_in = -ray.dir;

			// we don't need to check for `current_medium_id != -1` because in that case scatter is set to false.
			// do next event estimation for every scatter
			Spectrum nee_contrib = next_event_estimation(scene, ray.org, current_medium_id, vertex, dir_in, bounces,
														 rng, false, x, y);
			radiance += current_path_throughput * nee_contrib * sigma_s;

			Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
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

			dir_pdf = phase_func_pdf;
			nee_p_cache = ray.org;
			multi_trans_pdf = 1;
		} else {
			break;
		}

		Real rr_prob = 1;
		if (bounces >= scene.options.rr_depth) {
			rr_prob = fmin(max(current_path_throughput), Real(0.95));
			if (next_pcg32_real<Real>(rng) > rr_prob) {
				// Terminate the path
				break;
			}
		}
		bounces += 1;
		current_path_throughput /= rr_prob;
	}
	return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
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
	Real dir_pdf = 0; // in solid angle measure
	Spectrum nee_p_cache;
	Real multi_trans_pdf = 1;

	Real eta_scale = Real(1);

	int max_depth = scene.options.max_depth;
	if (debug(x, y)) {
		debug(x, y);
	}

	while(true) {
		bool scatter = false;

		std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
		Real t_hit = INF;
		PathVertex vertex;
		if (vertex_) {
			vertex = *vertex_;
			t_hit = distance(ray.org, vertex.position);
		}


		Real trans_pdf = 1;
		Real transmittance = 1;
		if (current_medium_id != -1) {
			Medium medium = scene.media[current_medium_id];
			Real sigma_a = get_sigma_a(medium, make_zero_spectrum()).x;
			Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;
			Real sigma_t = sigma_a + sigma_s;
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
			multi_trans_pdf *= trans_pdf;
		} else if (t_hit != INF) {
			ray.org = vertex.position;
		}
		current_path_throughput *= (transmittance/trans_pdf);

		// emission for phase function sampling
		// hit a surface
		if (!scatter && t_hit != INF  && is_light(scene.shapes[vertex.shape_id])) {

			// set emission if there is an intersection with a light source
			if (bounces == 0 ) {
				radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
			} else {
				// Need to account for next event estimation
				int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
				assert(light_id >= 0);
				const Light &light = scene.lights[light_id];
				PointAndNormal light_point{vertex.position, vertex.geometry_normal};

				Real pdf_nee = light_pmf(scene, light_id) *
							   pdf_point_on_light(light, light_point, nee_p_cache, scene);

				Real G = fmax(Real(0), -dot(ray.dir, vertex.geometry_normal)) /
						 distance_squared(nee_p_cache, vertex.position);
				Real pdf_phase = dir_pdf * multi_trans_pdf * G;
				Real weight = (pdf_phase * pdf_phase) / (pdf_phase * pdf_phase + pdf_nee * pdf_nee);
				radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * weight;
				if (debug(x, y)) {
					print(radiance, "NEE radiance in phase fn sampling in bounce " + std::to_string(bounces));
				}
			}
		}

		// index-matching interface, skip through it
		if (!scatter && t_hit != INF && vertex.material_id == -1) {
			current_medium_id = update_medium_id(ray, vertex, current_medium_id);
			ray = Ray{vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};

			bounces++;
			continue;
		}

		if (bounces >= max_depth - 1 and max_depth != -1) {
			break;
		}

		// radiance from nee will be added here
		// two cases arise
		// 1. Ray has been volumetrically scattered -> do NEE and phase function sampling
		// 2. Ray a hit a surface -> do NEE and BSDF sampling
		const Vector3 dir_in = -ray.dir;
		if (scatter && current_medium_id != -1) {
			Medium medium = scene.media[current_medium_id];
			Real sigma_s = get_sigma_s(medium, make_zero_spectrum()).x;
			// we don't need to check for `current_medium_id != -1` because in that case scatter is set to false.
			// do next event estimation for every scatter
			Spectrum nee_contrib = next_event_estimation(scene, ray.org, current_medium_id, vertex, dir_in, bounces, rng,
														 false, x, y);
			radiance += current_path_throughput * nee_contrib  * sigma_s;
			if (debug(x, y) ) {
				print(radiance, "NEE radiance in phase fn block in bounce " + std::to_string(bounces));
			}

			Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
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

			dir_pdf = phase_func_pdf;
			nee_p_cache = ray.org;
			multi_trans_pdf = 1;
		}
		else if (t_hit != INF) {
			// do NEE when we hit a surface
			// what if it is a light source? We don't care we still do NEE. The radiance from the light source will
			// dominate radiance of other bounces.
			// We don't need to check for `current_medium_id != -1` because in that case scatter is set to false.
			// do next event estimation for every scatter

			Spectrum nee_contrib = next_event_estimation(scene, ray.org, current_medium_id, vertex,
														 dir_in, bounces, rng, true, x, y);
			radiance += current_path_throughput * nee_contrib ;
			if (debug(x, y) ) {
				print(radiance, "NEE radiance in bsdf block in bounce " + std::to_string(bounces));
			}

			// bsdf sampling
			const Material &mat = scene.materials[vertex.material_id];

			Vector3 dir_view = -ray.dir;
			Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
			Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
			std::optional<BSDFSampleRecord> bsdf_sample_ =
				sample_bsdf(mat,
							dir_view,
							vertex,
							scene.texture_pool,
							bsdf_rnd_param_uv,
							bsdf_rnd_param_w);
			if (!bsdf_sample_) {
				break;
			}
			const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
			Vector3 next_dir = bsdf_sample.dir_out;
			// Update ray differentials & eta_scale
			if (bsdf_sample.eta == 0) {
				ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
			} else {
				ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
				eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
				current_medium_id = update_medium_id(ray, vertex, current_medium_id);
			}

			// Trace a ray towards bsdf_dir. Note that again we have
			// to have an "epsilon" tnear to prevent self intersection.
			ray = Ray{ray.org, next_dir, get_intersection_epsilon(scene), infinity<Real>()};

			Spectrum bsdf_eval = eval(mat, dir_view, next_dir, vertex, scene.texture_pool);
			Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, next_dir, vertex, scene.texture_pool);
			current_path_throughput *= (bsdf_eval/bsdf_pdf);

			dir_pdf = bsdf_pdf;
			nee_p_cache = ray.org;
			multi_trans_pdf = 1;
		}

		Real rr_prob = 1;
		if (bounces >= scene.options.rr_depth) {
			rr_prob = fmin(max((1 / eta_scale) * current_path_throughput), Real(0.95));
			if (next_pcg32_real<Real>(rng) > rr_prob) {
				// Terminate the path
				break;
			}
		}
		bounces += 1;
		current_path_throughput /= rr_prob;
	}
	return radiance;
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
