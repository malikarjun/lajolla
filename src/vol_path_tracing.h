#pragma once

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
		Ray vis_ray{p, dir_out,0, infinity<Real>()};

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
		Real trans_pdf = exp(-sigma_t * t);
		Real transmittance = exp(-sigma_t * t);
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

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
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
