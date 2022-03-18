#include "../microfacet.h"

Real get_Dm(Vector3 h, Vector2 alpha)
{
    return 1.0 / (c_PI * alpha.x * alpha.y * pow(pow((h.x / alpha.x), 2) + pow((h.y / alpha.y), 2) + h.z * h.z, 2));
}

Real get_G_metal(Vector3 omega, Vector2 alpha)
{
    double lambda = (sqrt(1 + (pow(omega.x * alpha.x, 2) + pow(omega.y * alpha.y, 2)) / (omega.z * omega.z)) - 1) / 2.0;
    return 1.0 / (1 + lambda);
}

Real get_Gm(Vector3 dir_in, Vector3 dir_out, Vector2 alpha)
{
    return get_G_metal(dir_in, alpha) * get_G_metal(dir_out, alpha);
}


Real rs(Real eta0, Real eta1, Real cosine0, Real cosine1)
{
    return (eta0 * cosine0 - eta1 * cosine1) / (eta0 * cosine0 + eta1 * cosine1);
}

Real rp(Real eta0, Real eta1, Real cosine0, Real cosine1)
{
    return (eta1 * cosine0 - eta0 * cosine1) / (eta0 * cosine1 + eta1 * cosine0);
}

Real ts(Real eta0, Real eta1, Real cosine0, Real cosine1)
{
    return (2 * eta0 * cosine0) / (eta0 * cosine0 + eta1 * cosine1);
}

Real tp(Real eta0, Real eta1, Real cosine0, Real cosine1)
{
    return (2 * eta0 * cosine0) / (eta0 * cosine1 + eta1 * cosine0);
}
//Reference: https://www.gamedev.net/tutorials/_/technical/graphics-programming-and-theory/thin-film-interference-for-computer-graphics-r2962/
std::pair<Real,Real> get_thinfilm_F(Vector3 dir_in, Vector3 n, float lambda, float thickness, Vector3 etas)
{
    lambda = lambda * (3.14 / 180.0);
    Real cosine0 = fabs(dot(dir_in, n));
    Real d10 = (etas[1] > etas[0]) ? 0 : c_PI;
    Real d12 = (etas[1] > etas[0]) ? 0 : c_PI;
    Real delta = d10 + d12;
    Real sine1 = pow(etas[0] / etas[1], 2) * (1 - pow(cosine0, 2));
    if (sine1 > 1)
    {
        return std::make_pair(1.0f,0.0f);
    }
    Real cosine1 = sqrt(1 - sine1);
    Real sine2 = pow(etas[0] / etas[2], 2) * (1 - pow(cosine0, 2));
    if (sine2 > 1)
    {
        return std::make_pair(1.0f, 0.0f);
    }
    Real cosine2 = sqrt(1 - sine2);
    Real alpha_s = rs(etas[1], etas[0], cosine1, cosine0) * rs(etas[1], etas[2], cosine1, cosine2);
    Real alpha_p = rp(etas[1], etas[0], cosine1, cosine0) * rp(etas[1], etas[2], cosine1, cosine2);
    Real beta_s = ts(etas[0], etas[1], cosine0, cosine1) * ts(etas[1], etas[2], cosine1, cosine2);
    Real beta_p = tp(etas[0], etas[1], cosine0, cosine1) * tp(etas[1], etas[2], cosine1, cosine2);
    Real phi = (2 * c_PI / lambda) * (2 * etas[1] * thickness * cosine1) + delta;
    Real ts = pow(beta_s, 2) / (pow(alpha_s, 2) - 2 * alpha_s * cos(phi) + 1);
    Real tp = pow(beta_p, 2) / (pow(alpha_p, 2) - 2 * alpha_p * cos(phi) + 1);
    Real beamRatio = (etas[2] * cosine2) / (etas[0] * cosine0);
    Real t = beamRatio * (ts + tp) / 2;
    return std::make_pair(1-t, t);
}


Vector2 alpha_aniso(double roughness, double aniso)
{
    double aspect = sqrt(1 - 0.9 * aniso);
    return Vector2(fmax(0.01, roughness * roughness / aspect), fmax(0.01, roughness * roughness * aspect));
}
Spectrum eval_op::operator()(const DisneyMetal& bsdf) const {
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
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real aniso = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }
    Vector2 alpha = alpha_aniso(roughness, aniso);

    //std::cout << dot(frame.n, dir_out)<<" "<< pow(dot(frame.n, dir_out), 5) << " " << (Real(1) - pow(dot(frame.n, dir_out), 5)) << std::endl;

    Vector3 etas{ 1.1,2.0,1.2 };
    Real thickness = 10;
    Real Fr = get_thinfilm_F(dir_in, frame.n, 650, thickness, etas).first;
    Real Fg = get_thinfilm_F(dir_in, frame.n, 510, thickness, etas).first;
    Real Fb = get_thinfilm_F(dir_in, frame.n, 475, thickness, etas).first;
    Spectrum Fm = Vector3{Fr,Fg,Fb} *(base_color + (Real(1) - base_color) * pow(Real(1) - dot(dir_out, half_vector), 5));
    Real Gm = get_Gm(to_local(frame, dir_in), to_local(frame, dir_out), alpha);
    Real Dm = get_Dm(to_local(frame, half_vector), alpha);




    return (Fm * Gm * Dm) / (4 * dot(frame.n, dir_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal& bsdf) const {
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

Vector3 sample_visible_normals_aniso(const Vector3& local_dir_in, Vector2 alpha, const Vector2& rnd_param) {
    // The incoming direction is in the "ellipsodial configuration" in Heitz's paper
    if (local_dir_in.z < 0) {
        // Ensure the input is on top of the surface.
        return -sample_visible_normals_aniso(-local_dir_in, alpha, rnd_param);
    }

    // Transform the incoming direction to the "hemisphere configuration".
    Vector3 hemi_dir_in = normalize(
        Vector3{ alpha.x * local_dir_in.x, alpha.y * local_dir_in.y, local_dir_in.z });

    // Parameterization of the projected area of a hemisphere.
    // First, sample a disk.
    Real r = sqrt(rnd_param.x);
    Real phi = 2 * c_PI * rnd_param.y;
    Real t1 = r * cos(phi);
    Real t2 = r * sin(phi);
    // Vertically scale the position of a sample to account for the projection.
    Real s = (1 + hemi_dir_in.z) / 2;
    t2 = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
    // Point in the disk space
    Vector3 disk_N{ t1, t2, sqrt(max(Real(0), 1 - t1 * t1 - t2 * t2)) };

    // Reprojection onto hemisphere -- we get our sampled normal in hemisphere space.
    Frame hemi_frame(hemi_dir_in);
    Vector3 hemi_N = to_world(hemi_frame, disk_N);

    // Transforming the normal back to the ellipsoid configuration
    return normalize(Vector3{ alpha.x * hemi_N.x, alpha.y * hemi_N.y, max(Real(0), hemi_N.z) });
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyMetal& bsdf) const {
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
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal& bsdf) const {
    return bsdf.base_color;
}
