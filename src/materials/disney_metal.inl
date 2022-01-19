#include "../microfacet.h"


Real sqr(Real val) {
	return val * val;
}


Real g_w(Vector3 w, Frame frame, Real alpha_x, Real alpha_y) {
	Vector3  w_l = to_local(frame, w);
	Real lambda_w = (sqrt(1 + (sqr(w_l.x * alpha_x) + sqr(w_l.y * alpha_y))/sqr(w_l.z) ) - 1)/2;
	return 1/(1 + lambda_w);
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

	// Clamp roughness to avoid numerical issues.
	roughness = std::clamp(roughness, Real(0.01), Real(1));

	Real h_dot_out = dot(half_vector, dir_out);
//	Real n_dot_out = dot(frame.n, dir_out);
	Real n_dot_in = dot(frame.n, dir_in);

	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = max(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = max(alpha_min, sqr(roughness) * aspect);

	Spectrum f_m = schlick_fresnel(base_color, h_dot_out);

	Real h_l_term_denom = sqr(sqr(h_l.x/ alpha_x) + sqr(h_l.y/ alpha_y) + sqr(h_l.z));
	Real d_m = 1/(c_PI * alpha_x * alpha_y * h_l_term_denom);

	Real g_m = g_w(dir_in, frame, alpha_x, alpha_y) * g_w(dir_out, frame, alpha_x, alpha_y);

	return f_m * d_m * g_m / (4 * n_dot_in);
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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


	Vector3 half_vector = normalize(dir_in + dir_out);
	// Flip half-vector if it's below surface
	if (dot(half_vector, frame.n) < 0) {
		half_vector = -half_vector;
	}

	Vector3 h_l = to_local(frame, half_vector);
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

	// Clamp roughness to avoid numerical issues.
	roughness = std::clamp(roughness, Real(0.01), Real(1));

	Real n_dot_in = dot(frame.n, dir_in);
	Real n_dot_out = dot(frame.n, dir_out);
	Real n_dot_h = dot(frame.n, half_vector);
	if (n_dot_out <= 0 || n_dot_h <= 0) {
		return 0;
	}

	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = max(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = max(alpha_min, sqr(roughness) * aspect);

	Real G = smith_masking_gtr2_anisotropic(to_local(frame, dir_in), alpha_x, alpha_y);

	// use anisotropic GTR2 instead of the isotropic one used in roughplastic
	Real h_l_term_denom = sqr(sqr(h_l.x/ alpha_x) + sqr(h_l.y/ alpha_y) + sqr(h_l.z));
	Real D = 1/(c_PI * alpha_x * alpha_y * h_l_term_denom);

//	Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness);
//	Real D = GTR2(n_dot_h, roughness);

	// (4 * cos_theta_v) is the Jacobian of the reflection
    return  (G * D) / (4 * n_dot_in);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

	// Sample from the specular lobe.
	// Convert the incoming direction to local coordinates
	Vector3 local_dir_in = to_local(frame, dir_in);
	Real roughness = eval(
		bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	// Clamp roughness to avoid numerical issues.
	roughness = std::clamp(roughness, Real(0.01), Real(1));
	Real alpha = roughness * roughness;
	Vector3 local_micro_normal =
		sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

	// Transform the micro normal to world space
	Vector3 half_vector = to_world(frame, local_micro_normal);
	// Reflect over the world space normal
	Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
	return BSDFSampleRecord{
		reflected,
		Real(0) /* eta */, roughness /* roughness */
	};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
