#include "../microfacet.h"


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

	// TODO : should h_dot_out be use abs
	Spectrum f_m = schlick_fresnel(base_color, h_dot_out);

	Real h_l_term_denom = sqr(sqr(h_l.x/ alpha_x) + sqr(h_l.y/ alpha_y) + sqr(h_l.z));
	Real d_m = 1/(c_PI * alpha_x * alpha_y * h_l_term_denom);

	Real G = smith_masking_gtr2_anisotropic(to_local(frame, dir_in), alpha_x, alpha_y) *
		smith_masking_gtr2_anisotropic(to_local(frame, dir_out), alpha_x, alpha_y);

	Spectrum val = f_m * d_m * G / (4 * fabs(n_dot_in));
	if (debug(x, y)) {
		print(val, "f_val");
	}
	return val;
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
//	roughness = std::clamp(roughness, Real(0.01), Real(1));

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
    Real val =  (G * D) / (4 * fabs(n_dot_in));

	if (false && debug(x, y)) {
		printf("Metal debugging!!\n");
		print(val);
		print(G, "G value");
		print(D, "D value");
		print(fabs(n_dot_in), "dot product ");
	}

	return val;
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
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

	// Clamp roughness to avoid numerical issues.
//	roughness = std::clamp(roughness, Real(0.01), Real(1));
	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = max(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = max(alpha_min, sqr(roughness) * aspect);
	Vector3 local_micro_normal =
		sample_visible_normals_anisotropic(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

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
