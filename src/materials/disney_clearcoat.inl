#include "../microfacet.h"

Real mask_shadow_c(Vector3 w, Frame frame, Real alpha = 0.25) {
	Vector3  w_l = to_local(frame, w);
	Real lambda_w = (sqrt(1 + (sqr(w_l.x * alpha) + sqr(w_l.y * alpha))/sqr(w_l.z) ) - 1)/2;
	return 1/(1 + lambda_w);
}

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
	Vector3 half_vector = normalize(dir_in + dir_out);
	// Flip half-vector if it's below surface
	if (dot(half_vector, frame.n) < 0) {
		half_vector = -half_vector;
	}
	Vector3 h_l = to_local(frame, half_vector);

	Real cleartcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
	cleartcoat_gloss = std::clamp(cleartcoat_gloss, Real(0.01), Real(1));

	Real h_dot_out = dot(half_vector, dir_out);
	Real n_dot_out = dot(frame.n, dir_out);
	Real n_dot_in = dot(frame.n, dir_in);
	Real n_dot_h = dot(frame.n, half_vector);

	if (n_dot_out <= 0 || n_dot_h <= 0) {
		return make_zero_spectrum();
	}

	Real eta = 1.5;
	Real alpha_g = (1 - cleartcoat_gloss) * 0.1 + cleartcoat_gloss * 0.001;
	Real alpha_g_2 = sqr(alpha_g);

	// TODO: should h_dot_out use abs
	Real f_c = schlick_fresnel(r_0(eta), h_dot_out);
	Real d_c = (alpha_g_2 - 1)/ (c_PI * log(alpha_g_2) * (1 + (alpha_g_2 - 1)* sqr(h_l.z)));
	Real g_c = mask_shadow_c(dir_in, frame) * mask_shadow_c(dir_out, frame);

	Real f_clearcoat = (f_c * d_c * g_c) / (4 * fabs(n_dot_in));

    return make_const_spectrum(1) * f_clearcoat;
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
	Vector3 half_vector = normalize(dir_in + dir_out);
	// Flip half-vector if it's below surface
	if (dot(half_vector, frame.n) < 0) {
		half_vector = -half_vector;
	}
	Vector3 h_l = to_local(frame, half_vector);

	Real cleartcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
	cleartcoat_gloss = std::clamp(cleartcoat_gloss, Real(0.01), Real(1));

	Real h_dot_out = dot(half_vector, dir_out);
	Real n_dot_out = dot(frame.n, dir_out);
//	Real n_dot_in = dot(frame.n, dir_in);
	Real n_dot_h = dot(frame.n, half_vector);

	if (n_dot_out <= 0 || n_dot_h <= 0) {
		return Real(0);
	}

	Real alpha_g = (1 - cleartcoat_gloss) * 0.1 + cleartcoat_gloss * 0.001;
	Real alpha_g_2 = sqr(alpha_g);

	Real d_c = (alpha_g_2 - 1)/ (c_PI * log(alpha_g_2) * (1 + (alpha_g_2 - 1)* sqr(h_l.z)));
    return d_c * n_dot_h / (4 * fabs(h_dot_out));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
	Real cleartcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
	cleartcoat_gloss = std::clamp(cleartcoat_gloss, Real(0.01), Real(1));
	Real alpha_g = (1 - cleartcoat_gloss) * 0.1 + cleartcoat_gloss * 0.001;
	Real alpha_g_2 = alpha_g * alpha_g;
	Real u0 = rnd_param_uv.x;
	Real u1 = rnd_param_uv.y;

	Real h_elevation = acos(sqrt((1 - pow(alpha_g_2, u0))/(1 - alpha_g_2)));
	Real h_azimuth = 2 * c_PI * u1;
	Vector3 h_l{fabs(sin(h_elevation)) * cos(h_azimuth), fabs(sin(h_elevation)) * sin(h_azimuth), cos(h_elevation)};

	// Transform the micro normal to world space
	Vector3 half_vector = to_world(frame, h_l);
	// Reflect over the world space normal
	Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
	return BSDFSampleRecord{
		reflected,
		Real(0) /* eta */, cleartcoat_gloss/* roughness */
	};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
