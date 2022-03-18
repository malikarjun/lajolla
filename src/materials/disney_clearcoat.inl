#include "../microfacet.h"

Real get_alpha(Real gloss)
{
    return (1 - gloss) * 0.1 + gloss * 0.001;
}

Real get_Dc(Real alpha, Vector3 hl)
{
    Real alpha_sq = alpha * alpha;
    return (alpha_sq - 1) / (c_PI * log(alpha_sq)*(1 + (alpha_sq - 1)*(hl.z*hl.z)));
}


Real get_R0(Real eta)
{
    return pow(eta - 1,2) / pow(eta + 1,2);
}

Real get_G_clearcoat(Vector3 omega)
{
    double lambda = (sqrt(1 + (pow(omega.x * 0.25, 2) + pow(omega.y * 0.25, 2)) / (omega.z * omega.z)) - 1) / 2.0;
    return 1.0 / (1 + lambda);
}

Real get_Gc(Vector3 dir_in, Vector3 dir_out)
{
    return get_G_clearcoat(dir_in) * get_G_clearcoat(dir_out);
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
    Vector3 h = normalize(dir_in + dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = get_alpha(gloss);
    Real Dc = get_Dc(alpha, to_local(frame, h));
    Real Fc = get_R0(1.5) + (1 - get_R0(1.5)) * pow((1 - abs(dot(h, dir_out))), 5);
    Real Gc = get_Gc(to_local(frame, dir_in), to_local(frame, dir_out));
    //std::cout << gloss << " " << alpha << " " << alpha << std::endl;
    //std::cout << Dc << " " << Fc << " " << Gc << std::endl;

    return make_zero_spectrum() + (Fc*Dc*Gc)/(4*abs(dot(frame.n,dir_in)));
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
    Vector3 h = normalize(dir_in + dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = get_alpha(gloss);
    Real Dc = get_Dc(alpha, to_local(frame, h));

    return Dc*abs(dot(frame.n,h))/(4*abs(dot(h,dir_out)));
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

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
