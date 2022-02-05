Real fresnel(Vector3 dir, Vector3 normal, Real f_90) {
  Real n_dot_dir = fabs(dot(dir, normal));
  return 1 + (f_90 - 1) * pow((1 - n_dot_dir), Real(5));
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // base diffuse
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
      half_vector = -half_vector;
    }
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

	// base diffuse
    Real h_dot_out = dot(half_vector, dir_out);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_in = dot(frame.n, dir_in);

    Real f_d90 = 0.5 + 2 * roughness * h_dot_out * h_dot_out;

    Real f_d_in = fresnel(dir_in, frame.n, f_d90);
    Real f_d_out = fresnel(dir_out, frame.n, f_d90);
    Spectrum f_base_diffuse = (base_color/c_PI) * f_d_in * f_d_out * fabs(n_dot_out);

    // subsurface
    Real f_ss90 = roughness * sqr(fabs(h_dot_out));
    Real f_ss_in = fresnel(dir_in, frame.n, f_ss90);
    Real f_ss_out = fresnel(dir_out, frame.n, f_ss90);
	Spectrum f_subsurface = ((1.25 * base_color )/ c_PI) *
								(f_ss_in * f_ss_out * (1/(fabs(n_dot_in) + fabs(n_dot_out)) - 0.5) + 0.5) * fabs(n_dot_out);

    return  (1 - subsurface) * f_base_diffuse + subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0 ||
            dot(vertex.geometry_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    // Homework 1: implement this!
	return BSDFSampleRecord{
		to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
		Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
