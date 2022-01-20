#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
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

	Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

	Real h_dot_out = dot(half_vector, dir_out);
	Real n_dot_out = dot(frame.n, dir_out);

	Real l_bc = luminance(base_color);
	Spectrum c_tint = make_const_spectrum(1);
	if (l_bc > 0) {
		c_tint = base_color/l_bc;
	}
	Spectrum c_sheen = (1 - sheen_tint) + sheen_tint * c_tint;
	Spectrum f_sheen = c_sheen * pow( (1 - fabs(h_dot_out)), 5) * fabs(n_dot_out);

    return f_sheen;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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
	return BSDFSampleRecord{
		to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
		Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
