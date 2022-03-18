#include "../microfacet.h"

Real get_Fg(Vector3 h, Vector3 dir_in, Vector3 dir_out, Real eta)
{
    //Real Rs = (dot(h, dir_in) - eta * dot(h, dir_out)) / (dot(h, dir_in) + eta * dot(h, dir_out));
    //Real Rp = (eta * dot(h, dir_in) - dot(h, dir_out)) / (eta * dot(h, dir_in) + dot(h, dir_out));
    //std::cout << dot(h, dir_in) << " "<< dot(h, dir_out) << std::endl;
    //return 0.5 * (Rs * Rs + Rp * Rp);
    Real h_dot_in = dot(h, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    return F;
}


Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector2 alpha = alpha_aniso(roughness, aniso);
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    // Homework 1: implement this!
    Real F = get_Fg(half_vector, dir_in, dir_out, eta);
    Real G = get_Gm(to_local(frame, dir_in), to_local(frame, dir_out), alpha);
    Real D = get_Dm(to_local(frame, half_vector), alpha);

    Vector3 etas{ 1.1,2.0,1.1 };
    Real thickness = 10;
    if (reflect)
    {
        Real Fr = get_thinfilm_F(dir_in, frame.n, 650, thickness, etas).first;
        Real Fg = get_thinfilm_F(dir_in, frame.n, 510, thickness, etas).first;
        Real Fb = get_thinfilm_F(dir_in, frame.n, 475, thickness, etas).first;
        return (base_color * F * (Vector3{ Fr,Fg,Fb }) * G * D) / (4 * abs(dot(frame.n, dir_in)));
    }
    else
    {
        Real Fr = get_thinfilm_F(dir_in, frame.n, 650, thickness, etas).second;
        Real Fg = get_thinfilm_F(dir_in, frame.n, 510, thickness, etas).second;
        Real Fb = get_thinfilm_F(dir_in, frame.n, 475, thickness, etas).second;
        Spectrum result = (sqrt(base_color) * (1 - F) * (Vector3 { Fr, Fg, Fb }) * D * G *
            abs(dot(half_vector, dir_out) * dot(half_vector, dir_in))) /
            (abs(dot(frame.n, dir_in)) * powf((dot(half_vector, dir_in) + eta * dot(half_vector, dir_out)), 2));
        //std::cout << F << std::endl;
        return result;
    }

    
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector2 alpha = alpha_aniso(roughness, aniso);
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = get_Fg(half_vector, dir_in, dir_out, eta);
    Real G = get_Gm(to_local(frame, dir_in), to_local(frame, dir_out), alpha);
    Real D = get_Dm(to_local(frame, half_vector), alpha);
    if (reflect) {
        return (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
    }
    else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometry_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector2 alpha = alpha_aniso(roughness, aniso);
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals_aniso(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{ reflected, Real(0) /* eta */, roughness };
    }
    else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
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
        Real h_dot_out = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{ refracted, eta, roughness };
    }

}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
