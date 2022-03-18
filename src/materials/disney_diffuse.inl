double F_D(Vector3 omega, double F_D90, const Frame& frame)
{
    return 1 + (F_D90 - 1) * pow((1 - dot(frame.n, omega)), 5);
}

double F_SS(Vector3 omega, double F_SS90, const Frame& frame)
{
    return 1 + (F_SS90 - 1) * pow((1 - dot(frame.n, omega)), 5);
}


Spectrum eval_op::operator()(const DisneyDiffuse& bsdf) const {
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
    double roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    double subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 h(normalize(dir_in+dir_out));
    double F_D90 = 0.5 + 2 * roughness * pow(dot(h, dir_out), 2);
    //

    Spectrum f_base = (base_color / c_PI)*F_D(dir_in, F_D90, frame)* F_D(dir_out, F_D90, frame)* fmax(dot(frame.n, dir_out), Real(0));

    double F_SS90 = roughness * pow(dot(h, dir_out), 2);

    Spectrum f_sub = 
        ((1.25 * base_color) / c_PI) * (F_SS(dir_in,F_SS90,frame)*
                                        F_SS(dir_out, F_SS90, frame)*
                                        (1/(fabs(dot(frame.n,dir_in))+fabs(dot(frame.n,dir_out))) - 0.5) + 0.5) * fmax(dot(frame.n, dir_out), Real(0));
    // Homework 1: implement this!
    return (1 - subsurface)* f_base + subsurface * f_sub;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse& bsdf) const {
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

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse& bsdf) const {
    if (dot(vertex.geometry_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    double roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(roughness) /* roughness */ };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse& bsdf) const {
    return bsdf.base_color;
}
