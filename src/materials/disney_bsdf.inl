#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

	Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);


	Real diffuseWeight = (1 - metallic) * (1 - specular_transmission);
	Real sheenWeight = (1 - metallic) * sheen;
	Real metalWeight = 1 - specular_transmission * (1 - metallic);
	Real clearcoatWeight = 0.25 * clearcoat;
	Real glassWeight = (1 - metallic) * specular_transmission;

//	sheenWeight = 0;
//	diffuseWeight = 0;
//	metalWeight = 0;
//	glassWeight = 1;
// 	clearcoatWeight = 0;


	Vector3 half_vector;

	if (reflect) {
		half_vector = normalize(dir_in + dir_out);
	} else {
		half_vector = normalize(dir_in + dir_out * eta);
	}
	if (dot(half_vector, frame.n) < 0) {
		half_vector = -half_vector;
	}
	Vector3 h_l = to_local(frame, half_vector);
	Real h_dot_out = dot(half_vector, dir_out);
	Real h_dot_in = dot(half_vector, dir_in);
	Real n_dot_out = dot(frame.n, dir_out);
	Real n_dot_in = dot(frame.n, dir_in);
	roughness = std::clamp(roughness, Real(0.01), Real(1));

	/* Diffuse  start */
	// base diffuse
	Real f_d90 = 0.5 + 2 * roughness * h_dot_out * h_dot_out;
	Real f_d_in = fresnel(dir_in, frame.n, f_d90);
	Real f_d_out = fresnel(dir_out, frame.n, f_d90);
	Spectrum f_base_diffuse = (base_color/c_PI) * f_d_in * f_d_out * fabs(n_dot_out);

	// subsurface
	Real f_ss90 = roughness * h_dot_out * h_dot_out;
	Real f_ss_in = fresnel(dir_in, frame.n, f_ss90);
	Real f_ss_out = fresnel(dir_out, frame.n, f_ss90);
	Spectrum f_subsurface = (1.25 * base_color / c_PI) *
							(f_ss_in * f_ss_out * (1/(fabs(n_dot_in) + fabs(n_dot_out)) - 0.5) + 0.5) * fabs(n_dot_out);
	Spectrum f_diffuse =  (1 - subsurface) * f_base_diffuse + subsurface * f_subsurface;
	/* Diffuse  end */

	/* Sheen start */
	Real l_bc = luminance(base_color);
	Spectrum c_tint = make_const_spectrum(1);
	if (l_bc > 0) {
		c_tint = base_color/l_bc;
	}
	Spectrum c_sheen = (1 - sheen_tint) + sheen_tint * c_tint;
	Spectrum f_sheen = c_sheen * pow( (1 - fabs(h_dot_out)), 5) * fabs(n_dot_out);
	/* Sheen end */

	/* Metal start */
	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = max(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = max(alpha_min, sqr(roughness) * aspect);

	// TODO : should h_dot_out be use abs
	Spectrum k_s = (1 - specular_tint) + specular_tint * c_tint;
	Spectrum c_0 = specular * r_0(eta) * (1 - metallic) * k_s + metallic * base_color;
	Spectrum f_m = schlick_fresnel(c_0, abs(h_dot_out));

	Real h_l_term_denom = sqr(sqr(h_l.x/ alpha_x) + sqr(h_l.y/ alpha_y) + sqr(h_l.z));
	Real D = 1/(c_PI * alpha_x * alpha_y * h_l_term_denom);

	Real G = smith_masking_gtr2_anisotropic(to_local(frame, dir_in), alpha_x, alpha_y) *
			 smith_masking_gtr2_anisotropic(to_local(frame, dir_out), alpha_x, alpha_y);

	Spectrum f_metal = f_m * D * G / (4 * fabs(n_dot_in));
	/* Metal end */

	/* Glass start */
	Real F = fresnel_dielectric(h_dot_in, eta);
	Spectrum f_glass;
	if (reflect) {
		f_glass = base_color * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
	} else {
		f_glass = sqrt(base_color) * ( (1 - F) * D * G * fabs(h_dot_out * h_dot_in)) /
			   (fabs(n_dot_in) * sqr(h_dot_in + eta * h_dot_out));
	}
	/* Glass end */

	/* Clearcoat start */
	Real eta_cc = 1.5;
	Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
	Real alpha_g_2 = sqr(alpha_g);

	Real f_c = schlick_fresnel(r_0(eta_cc), h_dot_out);
	Real d_c = (alpha_g_2 - 1)/ (c_PI * log(alpha_g_2) * (1 + (alpha_g_2 - 1)* sqr(h_l.z)));
	Real g_c = mask_shadow_c(dir_in, frame) * mask_shadow_c(dir_out, frame);

	Spectrum f_clearcoat = make_const_spectrum(1) * (f_c * d_c * g_c) / (4 * fabs(n_dot_in));
	/* Clearcoat end */

	if (dot(frame.n, dir_in) < 0 || dot(frame.n, dir_out) < 0) {
		diffuseWeight = 0;
		sheenWeight = 0;
		metalWeight = 0;
		clearcoatWeight = 0;
	}

//	return f_diffuse;
// 	return f_metal;
//	return f_clearcoat; TODO: not working
//	return f_glass;
//	return f_sheen;

	if (reflect) {
		return f_diffuse * diffuseWeight + f_sheen * sheenWeight + f_metal * metalWeight +
			   f_glass * glassWeight + f_clearcoat * clearcoatWeight;
	} else {
		return  f_glass * glassWeight;
	}
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
	Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

	Real diffuseWeight = (1 - metallic) * (1 - specular_transmission);
	Real sheenWeight = (1 - metallic) * sheen;
	Real metalWeight = 1 - specular_transmission * (1 - metallic);
	Real clearcoatWeight = 0.25 * clearcoat;
	Real glassWeight = (1 - metallic) * specular_transmission;


//	diffuseWeight = 0;
//	metalWeight = 0;
//	glassWeight = 1;
// 	clearcoatWeight = 0;


	Real totalWeight = diffuseWeight + metalWeight + clearcoatWeight + glassWeight;
	diffuseWeight /= totalWeight;
	metalWeight /= totalWeight;
	clearcoatWeight /= totalWeight;
	glassWeight /= totalWeight;

	Vector3 half_vector;

	if (reflect) {
		half_vector = normalize(dir_in + dir_out);
	} else {
		half_vector = normalize(dir_in + dir_out * eta);
	}
	if (dot(half_vector, frame.n) < 0) {
		half_vector = -half_vector;
	}
	Vector3 h_l = to_local(frame, half_vector);
	Real h_dot_out = dot(half_vector, dir_out);
	Real h_dot_in = dot(half_vector, dir_in);
	Real n_dot_out = dot(frame.n, dir_out);
	Real n_dot_in = dot(frame.n, dir_in);
	Real n_dot_h = dot(frame.n, half_vector);
	roughness = std::clamp(roughness, Real(0.01), Real(1));

	/* Diffuse start */
	Real pdf_diffuse = fmax(n_dot_out, Real(0)) / c_PI;
	/* Diffuse  end */

	/* Metal start */
	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = max(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = max(alpha_min, sqr(roughness) * aspect);

	Real G = smith_masking_gtr2_anisotropic(to_local(frame, dir_in), alpha_x, alpha_y);

	Real h_l_term_denom = sqr(sqr(h_l.x/ alpha_x) + sqr(h_l.y/ alpha_y) + sqr(h_l.z));
	Real D = 1/(c_PI * alpha_x * alpha_y * h_l_term_denom);

	Real pdf_metal = (G * D) / (4 * fabs(n_dot_in));
	/* Metal  end */

	/* Glass start */
	Real F = fresnel_dielectric(h_dot_in, eta);
	Real pdf_glass;
	if (reflect) {
		pdf_glass =  (F * D * G) / (4 * fabs(n_dot_in));
	} else {
		Real sqrt_denom = h_dot_in + eta * h_dot_out;
		Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
		pdf_glass =  (1 - F) * D * G * fabs(dh_dout * h_dot_in / n_dot_in);
	}
	/* Glass end */

	/* Clearcoat start */
	Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
	Real alpha_g_2 = sqr(alpha_g);

	Real d_c = (alpha_g_2 - 1)/ (c_PI * log(alpha_g_2) * (1 + (alpha_g_2 - 1)* sqr(h_l.z)));
	Real pdf_clearcoat = d_c * n_dot_h / (4 * fabs(h_dot_out));
	/* Clearcoat end */

	if (dot(frame.n, dir_in) < 0 || dot(frame.n, dir_out) < 0) {
		// Our incoming ray is coming from inside,
		// so the probability of sampling the glass lobe is 1 if glass_prob is not 0.
		diffuseWeight = 0;
		metalWeight = 0;
		clearcoatWeight = 0;
		if (glassWeight > 0) {
			glassWeight = 1;
		}
	}

//	return pdf_diffuse;
//	return pdf_clearcoat;
//	return pdf_glass;
	if (reflect) {
		return pdf_diffuse * diffuseWeight + pdf_metal * metalWeight + pdf_glass * glassWeight +
		pdf_clearcoat * clearcoatWeight;
	} else {
		return  pdf_glass * glassWeight;
	}
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
	Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);


	Real diffuseWeight = (1 - metallic) * (1 - specular_transmission);
	Real metalWeight = 1 - specular_transmission * (1 - metallic);
	Real clearcoatWeight = 0.25 * clearcoat;
	Real glassWeight = (1 - metallic) * specular_transmission;

//	diffuseWeight = 0;
//	metalWeight = 0;
//	glassWeight = 1;
// 	clearcoatWeight = 0;

	Real n_dot_in = dot(frame.n, dir_in);
	roughness = std::clamp(roughness, Real(0.01), Real(1));

	Real aspect = sqrt(1 - 0.9 * anisotropic);
	Real alpha_min = 0.0001;
	Real alpha_x = fmax(alpha_min, sqr(roughness)/ aspect);
	Real alpha_y = fmax(alpha_min, sqr(roughness) * aspect);

	Vector3 local_dir_in = to_local(frame, dir_in);
	Vector3 local_micro_normal =
		sample_visible_normals_anisotropic(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

	Real totalWeight = diffuseWeight  + metalWeight + glassWeight + clearcoatWeight;
	diffuseWeight /= totalWeight;
	metalWeight /= totalWeight;
	glassWeight /= totalWeight;
	clearcoatWeight /= totalWeight;

	Real cdfDiffuse = diffuseWeight;
	Real cdfMetal = cdfDiffuse + metalWeight;
	Real cdfGlass = cdfMetal + glassWeight;
//	Real cdfClearcoat = cdfGlass + clearcoatWeight; // this should be 1


	if (dot(frame.n, dir_in) >= 0) {
		// choose between diffuse, metallic, glass, and clearcoat
		if (rnd_param_w <= cdfDiffuse) {
			if (debug(x, y))
				printf("sampling diffuse!!!\n");
			return BSDFSampleRecord{
				to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
				Real(0), roughness };
		}
		else if (rnd_param_w <= cdfMetal) {
			// Transform the micro normal to world space
			if (debug(x, y))
				printf("sampling metal at (%d, %d)!!!\n", x, y);
			Vector3 half_vector = to_world(frame, local_micro_normal);
			// Reflect over the world space normal
			Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
			return BSDFSampleRecord{ reflected,Real(0) , roughness};
		}
		else if (rnd_param_w <= cdfGlass) {
			// choose between reflection & refraction
			Vector3 half_vector = to_world(frame, local_micro_normal);
			// Flip half-vector if it's below surface
			if (dot(half_vector, frame.n) < 0) {
				half_vector = -half_vector;
			}
			Real h_dot_in = dot(half_vector, dir_in);
			Real F = fresnel_dielectric(h_dot_in, eta);

//			printf("Above surface Fresnel is %f \n", F);
			const Real &scaled_w = (rnd_param_w - cdfMetal) / (cdfGlass - cdfMetal);
			if (scaled_w <= F) {
				nreflect += 1;
				if (debug(x, y))
					printf("sampling glass reflection at (%d, %d)!!!\n", x, y);
				// Reflection
				Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
				// set eta to 0 since we are not transmitting
				return BSDFSampleRecord{reflected, Real(0), roughness};
			} else {
				// Refraction
				nrefract += 1;
				if (debug(x, y))
					printf("sampling glass refraction at (%d, %d)!!!\n", x, y);
				Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
				if (h_dot_out_sq <= 0) {
					// Total internal reflection
					// This shouldn't really happen, as F will be 1 in this case.
					return {};
				}
				// flip half_vector if needed
				if (h_dot_in < 0) {
					half_vector = -half_vector;
				}
				Real h_dot_out= sqrt(h_dot_out_sq);
				Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
				return BSDFSampleRecord{refracted, eta, roughness};
			}
		}
		else {
			if (debug(x, y))
				printf("sampling clearcoat!!!\n");
			Real alpha_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
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
			return BSDFSampleRecord{reflected,
									Real(0) , roughness};
//			return BSDFSampleRecord{reflected,
//									Real(0) , fmax(0, sqrt(alpha_g))};
		}
	}
	else
	{
		// sample glass
		// choose between reflection & refraction
		Vector3 half_vector = to_world(frame, local_micro_normal);
		// Flip half-vector if it's below surface
		if (dot(half_vector, frame.n) < 0) {
			half_vector = -half_vector;
		}
		Real h_dot_in = dot(half_vector, dir_in);
		Real F = fresnel_dielectric(h_dot_in, eta);

//		printf("Below surface Fresnel is %f \n", F);
		if (rnd_param_w <= F) {
			nreflect += 1;
			if (debug(x, y))
				printf("sampling glass reflection at (%d, %d)!!!\n", x, y);
			// Reflection
			Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
			// set eta to 0 since we are not transmitting
			return BSDFSampleRecord{reflected, Real(0) , roughness};
		}

		else {
			nrefract += 1;
			if (debug(x, y))
				printf("sampling glass refraction!!!\n");
			// Refraction
			Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
			if (h_dot_out_sq <= 0) {
				// Total internal reflection
				// This shouldn't really happen, as F will be 1 in this case.
				return {};
			}
			// flip half_vector if needed
			if (h_dot_in < 0) {
				half_vector = -half_vector;
			}
			Real h_dot_out= sqrt(h_dot_out_sq);
			Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
			return BSDFSampleRecord{refracted, eta, roughness};
		}
	}
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
