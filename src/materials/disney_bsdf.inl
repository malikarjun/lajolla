#include "../microfacet.h"

//-------------------------diffuse------------------------------//
Real diffuse_pdf(
    const Frame& frame, 
    const Vector3& dir_in, 
    const Vector3& dir_out, 
    const DisneyBSDF& bsdf,
    const PathVertex& vertex,
    const TexturePool& texture_pool
    ) {
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}
std::optional<BSDFSampleRecord> diffuse_sample(
    const Frame& frame,
    const Vector3& dir_in,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf,
    const Vector2& rnd_param_uv)
{
    if (dot(frame.n, dir_in) < 0) {
        // No light below the surface
        return {};
    }

    double roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(roughness) /* roughness */ };
}

Spectrum diffuse_eval(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf)
{
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 h(normalize(dir_in + dir_out));
    Real F_D90 = 0.5 + 2 * roughness * pow(dot(h, dir_out), 2);
    //

    Spectrum f_base = (base_color / c_PI) * F_D(dir_in, F_D90, frame) * F_D(dir_out, F_D90, frame) * fmax(dot(frame.n, dir_out), Real(0));

    Real F_SS90 = roughness * pow(dot(h, dir_out), 2);

    Spectrum f_sub =
        ((1.25 * base_color) / c_PI) * (F_SS(dir_in, F_SS90, frame) *
            F_SS(dir_out, F_SS90, frame) *
            (1 / (fabs(dot(frame.n, dir_in)) + fabs(dot(frame.n, dir_out))) - 0.5) + 0.5) * fmax(dot(frame.n, dir_out), Real(0));
    Spectrum result = (1 - subsurface) * f_base + subsurface * f_sub;
    //if (f_sub.x != f_sub.x)
    //{
    //    return make_zero_spectrum();
    //}
    return result;
}
//-------------------------end of diffuse------------------------------//

//-------------------------sheen------------------------------//
Spectrum sheen_eval(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf)
{
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Homework 1: implement this!
    Spectrum Ctint = luminance(base_color) > 0 ? base_color / luminance(base_color) : Vector3(1.0, 1.0, 1.0);
    Spectrum Csheen = (1 - sheen_tint) + sheen_tint * Ctint;
    return  Csheen * pow(1 - fabs(dot(half_vector, dir_out)), 5) * fabs(dot(frame.n, dir_out));
}
//-------------------------end of sheen------------------------------//
//-------------------------clearcoat------------------------------//
Real clearcoat_pdf(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const DisneyBSDF& bsdf,
    const PathVertex& vertex,
    const TexturePool& texture_pool
) {
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    Vector3 h = normalize(dir_in + dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = get_alpha(gloss);
    Real Dc = get_Dc(alpha, to_local(frame, h));
    //std::cout << abs(dot(frame.n, h)) << "\n";
    return  (Dc * abs(dot(frame.n, h))) / (4 * abs(dot(h, dir_out)));
}
std::optional<BSDFSampleRecord> clearcoat_sample(
    const Frame& frame,
    const Vector3& dir_in,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf,
    const Vector2& rnd_param_uv)
{
    if (dot(frame.n, dir_in) < 0) {
        // No light below the surface
        return {};
    }

    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = get_alpha(gloss);
    Real h_ele = acos(sqrt((1 - pow((alpha * alpha), 1 - rnd_param_uv.x)) / (1 - alpha * alpha)));
    Real h_azi = 2 * c_PI * rnd_param_uv.y;

    Vector3 hl(sin(h_ele) * cos(h_azi), sin(h_ele) * sin(h_azi), cos(h_ele));

    Vector3 half_vector = to_world(frame, hl);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);

    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, alpha /* roughness */ };
}

Spectrum clearcoat_eval(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf)
{
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }

    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = get_alpha(gloss);
    Real Dc = get_Dc(alpha, to_local(frame, h));
    Real Fc = get_R0(1.5) + (1 - get_R0(1.5)) * pow((1 - abs(dot(h, dir_out))), 5);
    Real Gc = get_Gc(to_local(frame, dir_in), to_local(frame, dir_out));
    //std::cout << gloss << " " << alpha << " " << alpha << std::endl;
    //std::cout << Dc << " " << Fc << " " << Gc << std::endl;

    return make_zero_spectrum() + (Fc * Dc * Gc) / (4 * abs(dot(frame.n, dir_in)));
}
//-------------------------end of clearcoat------------------------------//

//-------------------------metal------------------------------//



Real metal_pdf(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const DisneyBSDF& bsdf,
    const PathVertex& vertex,
    const TexturePool& texture_pool
) {
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return 0;
    }

    // Homework 1: implement this!
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector2 alpha = alpha_aniso(roughness, aniso);
    Real G = get_G_metal(to_local(frame, dir_in), alpha);
    Real D = get_Dm(to_local(frame, half_vector), alpha);

    return (G * D) / (4 * n_dot_in);
}
std::optional<BSDFSampleRecord> metal_sample(
    const Frame& frame,
    const Vector3& dir_in,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf,
    const Vector2& rnd_param_uv)
{
    if (dot(frame.n, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);
    // Homework 1: implement this!
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector2 alpha = alpha_aniso(roughness, aniso);
    Vector3 local_micro_normal =
        sample_visible_normals_aniso(local_dir_in, alpha, rnd_param_uv);

    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(eta) /* eta */, roughness /* roughness */
    };
}

Spectrum metal_eval(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf)
{
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);
    // Homework 1: implement this!
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }
    Vector2 alpha = alpha_aniso(roughness, aniso);

    //std::cout << dot(frame.n, dir_out)<<" "<< pow(dot(frame.n, dir_out), 5) << " " << (Real(1) - pow(dot(frame.n, dir_out), 5)) << std::endl;
    Spectrum Ctint = luminance(base_color) > 0 ? base_color / luminance(base_color) : Vector3(1.0, 1.0, 1.0);
    Spectrum Ks = (1 - specular_tint) + specular_tint * Ctint;
    Spectrum C0 = specular * (get_R0(eta)) * (1 - metallic) * Ks + metallic * base_color;

    //Spectrum Fm = C0 + (Real(1) - C0) * pow(Real(1) - fabs(dot(dir_out, half_vector)), 5);
    Spectrum Fm = C0 + (1.0 - C0) * pow(1.0 - dot(dir_out, half_vector), 5);
    //std::cout << C0 << "//" << Ctint << "//" << Ctint << "//" << base_color << std::endl;
    Real Gm = get_Gm(to_local(frame, dir_in), to_local(frame, dir_out), alpha);
    Real Dm = get_Dm(to_local(frame, half_vector), alpha);

    return (Fm * Gm * Dm) / (4 * dot(frame.n, dir_in));
}
//-------------------------end of metal------------------------------//

//-------------------------glass------------------------------//
Real glass_pdf(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const DisneyBSDF& bsdf,
    const PathVertex& vertex,
    const TexturePool& texture_pool
) {
    bool reflect = dot(frame.n, dir_in) *
        dot(frame.n, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    // Homework 1: implement this!
    Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
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
std::optional<BSDFSampleRecord> glass_sample(
    const Frame& frame,
    const Vector3& dir_in,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf,
    const Vector2& rnd_param_uv,
    const Real& rnd_param_w)
{
    // Homework 1: implement this!
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(frame.n, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
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

Spectrum glass_eval(
    const Frame& frame,
    const Vector3& dir_in,
    const Vector3& dir_out,
    const PathVertex& vertex,
    const TexturePool& texture_pool,
    const DisneyBSDF& bsdf)
{
    bool reflect = dot(frame.n, dir_in) *
        dot(frame.n, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
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
    if (reflect)
    {
        return (base_color * F * G * D) / (4 * abs(dot(frame.n, dir_in)));
    }
    else
    {
        Spectrum result = (sqrt(base_color) * (1 - F) * D * G *
            abs(dot(half_vector, dir_out) * dot(half_vector, dir_in))) /
            (abs(dot(frame.n, dir_in)) * powf((dot(half_vector, dir_in) + eta * dot(half_vector, dir_out)), 2));
        //std::cout << F << std::endl;
        return result;
    }
}
//-------------------------end of glass------------------------------//

//******************************************************************************************************************//
Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometry_normal, dir_in) *
                   dot(vertex.geometry_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; 
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_sheen = (1 - metallic) * sheen;
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    //std::cout << w_glass << "\n";
    //Real sw_sum = w_diffuse + w_sheen + w_metal + w_clearcoat + w_glass;
    //w_diffuse /= sw_sum;
    //w_sheen /= sw_sum;
    //w_metal /= sw_sum;
    //w_clearcoat /= sw_sum;
    //w_glass /= sw_sum;

    Spectrum result;
    Spectrum glass_result;
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0)
    {
        result = w_glass * glass_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf);
    }
    else
    {
        if (reflect)
        {
            result =
                w_diffuse * diffuse_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
                w_clearcoat * clearcoat_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
                w_sheen * sheen_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
                w_metal * metal_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
                w_glass * glass_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf);
        }
        else
        {
            result = w_glass * glass_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf);
        }
    }
    //result = w_metal * metal_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
    //    w_clearcoat * clearcoat_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
    //    w_diffuse * diffuse_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf) +
    //    w_sheen * sheen_eval(frame, dir_in, dir_out, vertex, texture_pool, bsdf);

    return result;
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

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sw_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real sw_clearcoat = 0.25 * clearcoat;
    Real sw_metal = (1 - specular_transmission * (1 - metallic));
    Real sw_glass = (1 - metallic) * specular_transmission;
    Real sw_sum = sw_diffuse + sw_clearcoat + sw_metal +sw_glass;
    sw_diffuse /= sw_sum;
    sw_clearcoat /= sw_sum;
    sw_metal /= sw_sum;
    sw_glass /= sw_sum;
    if (dot(frame.n, dir_in) < 0 ||
        dot(frame.n, dir_out) < 0)
    {
        return glass_pdf(frame, dir_in, dir_out, bsdf, vertex, texture_pool);
    }
    return  sw_clearcoat* clearcoat_pdf(frame, dir_in, dir_out, bsdf, vertex, texture_pool)+
            sw_metal * metal_pdf(frame, dir_in, dir_out, bsdf, vertex, texture_pool) + 
            sw_diffuse * diffuse_pdf(frame, dir_in, dir_out, bsdf, vertex, texture_pool) +
            sw_glass * glass_pdf(frame, dir_in, dir_out, bsdf, vertex, texture_pool);

}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometry_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    
    Real sw_clearcoat = 0.25 * clearcoat;
    Real sw_metal = (1 - specular_transmission * (1 - metallic));
    Real sw_glass = (1 - metallic) * specular_transmission;
    Real sw_diffuse =(1 - specular_transmission) * (1 - metallic);//////////////
    Real sw_sum = sw_diffuse + sw_clearcoat + sw_metal +sw_glass;//only for testing because adding compoent step by step
    sw_diffuse /= sw_sum;
    sw_clearcoat /= sw_sum;
    sw_metal /= sw_sum;
    sw_glass /= sw_sum;

    Real cdf_clearcoat = sw_clearcoat;
    Real cdf_metal = cdf_clearcoat + sw_metal;
    Real cdf_glass = cdf_metal+ sw_glass;

    Real rnd = (std::rand() / Real(RAND_MAX));
    //std::cout << 0.25 * clearcoat << " " << (1 - specular_transmission * (1 - metallic)) << " " << (1 - metallic) * specular_transmission <<" "<< (1 - specular_transmission) * (1 - metallic) << std::endl;
    if (dot(dir_in, frame.n) <= 0)
    {
        return glass_sample(frame, dir_in, vertex, texture_pool, bsdf, rnd_param_uv, rnd_param_w);
    }
    if (rnd  < cdf_clearcoat)
    {
        return clearcoat_sample(frame, dir_in, vertex, texture_pool, bsdf, rnd_param_uv);
    }
    else if (rnd < cdf_metal)
    {
        return metal_sample(frame, dir_in, vertex, texture_pool, bsdf, rnd_param_uv);
    }
    else if (rnd < cdf_glass)
    {
        return glass_sample(frame, dir_in, vertex, texture_pool, bsdf, rnd_param_uv, rnd_param_w);
    }
    else
    {
        return diffuse_sample(frame, dir_in, vertex, texture_pool, bsdf, rnd_param_uv);
    }

}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
