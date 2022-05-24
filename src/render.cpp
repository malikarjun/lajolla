#include "render.h"
#include "intersection.h"
#include "material.h"
#include "parallel.h"
#include "path_tracing.h"
#include "vol_path_tracing.h"
#include "photon_mapping.h"
#include "pcg.h"
#include "progress_reporter.h"
#include "scene.h"
#include "transform.h"


/// Render auxiliary buffers e.g., depth.
Image3 aux_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    Image4 img4(w, h);
/*    for (int k = 0; k < 4; ++k) {
        for (int x = 0; x < img.width; ++x) {
            for (int y = 0; y < img.height; ++y) {
                image4(x, y) = Vector4{1, 2, 3, 4};
            }
        }
    }

    save_npy(image4, "tmp.npy");*/

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Ray ray = sample_primary(scene.camera, Vector2((x + Real(0.5)) / w, (y + Real(0.5)) / h));
                RayDifferential ray_diff = init_ray_differential(w, h);
                if (std::optional<PathVertex> vertex = intersect(scene, ray, ray_diff)) {
                    // Real dist = distance(vertex->position, ray.org);
                    Vector3 color{0, 0, 0};
                    Vector4 color4 {-1, -1, -1, -1};
                    if (scene.options.integrator == Integrator::Depth) {
                        Real z_eye_space = xform_point(scene.camera.world_to_cam, vertex->position)[2];
                        color = Vector3{z_eye_space, z_eye_space, z_eye_space};

                        color4.x = vertex->shape_id;
                        color4.y = vertex->primitive_id;
                        color4.z = vertex->st.x;
                        color4.w = vertex->st.y;

                        // color = Vector3{dist, dist, dist};
                    } else if (scene.options.integrator == Integrator::ShadingNormal) {
                        // color = (vertex->shading_frame.n + Vector3{1, 1, 1}) / Real(2);
                        color = vertex->shading_frame.n;
                    } else if (scene.options.integrator == Integrator::MeanCurvature) {
                        Real kappa = vertex->mean_curvature;
                        color = Vector3{kappa, kappa, kappa};
                    } else if (scene.options.integrator == Integrator::RayDifferential) {
                        color = Vector3{ray_diff.radius, ray_diff.spread, Real(0)};
                    } else if (scene.options.integrator == Integrator::MipmapLevel) {
                        const Material &mat = scene.materials[vertex->material_id];
                        const TextureSpectrum &texture = get_texture(mat);
                        auto *t = std::get_if<ImageTexture<Spectrum>>(&texture);
                        if (t != nullptr) {
                            const Mipmap3 &mipmap = get_img3(scene.texture_pool, t->texture_id);
                            Vector2 uv{modulo(vertex->uv[0] * t->uscale, Real(1)),
                                       modulo(vertex->uv[1] * t->vscale, Real(1))};
                            // ray_diff.radius stores approximatedly dpdx,
                            // but we want dudx -- we get it through
                            // dpdx / dpdu
                            Real footprint = vertex->uv_screen_size;
                            Real scaled_footprint = max(get_width(mipmap), get_height(mipmap)) *
                                                    max(t->uscale, t->vscale) * footprint;
                            Real level = log2(max(scaled_footprint, Real(1e-8f)));
                            color = Vector3{level, level, level};
                        }
                    }
                    img(x, y) = color;
                    img4(x, y) = color4;
                } else {
                    img(x, y) = Vector3{0, 0, 0};
                }
            }
        }
    }, Vector2i(num_tiles_x, num_tiles_y));

    if (scene.options.integrator == Integrator::Depth) {
        // convert eye space/world space depth into clip space depth needed for A-SVGF
        // TODO: not sure if near clip plane is 1
        // FIXME: Ans. Checked by debugging near clip plane is 1, but it is not enforced due to lack of camera frustum.
        Real n = 1, f = -1;
        //
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                f = std::max(f, img(x, y)[0]);
            }
        }
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                img(x, y) = (img(x, y) + Real(1))/(f - n);
            }
        }

        std::string base_path = replace(scene.output_filename, "_depth.exr", "");
        save_mat(scene.camera.cam_to_sample * scene.camera.world_to_cam, base_path + "_viewproj.npy");
        save_npy(img4, base_path + "_vbuffer.npy");
        save_npy(scene.model_mats, base_path + "_model_mats.npy");
        save_txt(scene.model_fnames, base_path + "_model_fnames.txt");
    }

    return img;
}

Image3 path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;
    #undef DEBUG
    #ifdef DEBUG

		ProgressReporter reporter(h*w);
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				Spectrum radiance = make_zero_spectrum();
				int spp = scene.options.samples_per_pixel;
				pcg32_state rng = init_pcg32(x*w + y);
                Real depth = 0;
                Vector3 normal;

//                if (debug(x, y)) {
//                    debug(0, 0);
//                }

				for (int s = 0; s < spp; s++) {
					Spectrum L = path_tracing(scene, x, y, rng, depth, normal);
					if (isfinite(L)) {
						// Hacky: exclude NaNs in the rendering.
						radiance += L;
					}
				}
				img(x, y) = radiance / Real(spp);
/*                depth_buffer[x][y] = depth;
                normal_buffer[x][y] = normal;*/
//				reporter.update(1);
			}

		}

		reporter.done();
    #else
		ProgressReporter reporter(num_tiles_x * num_tiles_y);
		parallel_for([&](const Vector2i &tile) {
			// Use a different rng stream for each thread.
			pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
			int x0 = tile[0] * tile_size;
			int x1 = min(x0 + tile_size, w);
			int y0 = tile[1] * tile_size;
			int y1 = min(y0 + tile_size, h);
			for (int y = y0; y < y1; y++) {
				for (int x = x0; x < x1; x++) {
					Spectrum radiance = make_zero_spectrum();
					int spp = scene.options.samples_per_pixel;
					for (int s = 0; s < spp; s++) {
						radiance += path_tracing(scene, x, y, rng);
					}
					img(x, y) = radiance / Real(spp);
				}
			}
			reporter.update(1);
		}, Vector2i(num_tiles_x, num_tiles_y));
		reporter.done();
    #endif



    return img;
}

Image3 vol_path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    auto f = vol_path_tracing;
    if (scene.options.vol_path_version == 1) {
        f = vol_path_tracing_1;
    } else if (scene.options.vol_path_version == 2) {
        f = vol_path_tracing_2;
    } else if (scene.options.vol_path_version == 3) {
        f = vol_path_tracing_3;
    } else if (scene.options.vol_path_version == 4) {
        f = vol_path_tracing_4;
    } else if (scene.options.vol_path_version == 5) {
        f = vol_path_tracing_5;
    } else if (scene.options.vol_path_version == 6) {
        f = vol_path_tracing;
    }

	if (1) {
		ProgressReporter reporter(num_tiles_x * num_tiles_y);
		parallel_for([&](const Vector2i &tile) {
			// Use a different rng stream for each thread.
			pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
			int x0 = tile[0] * tile_size;
			int x1 = min(x0 + tile_size, w);
			int y0 = tile[1] * tile_size;
			int y1 = min(y0 + tile_size, h);
			for (int y = y0; y < y1; y++) {
				for (int x = x0; x < x1; x++) {
					Spectrum radiance = make_zero_spectrum();
					int spp = scene.options.samples_per_pixel;
					for (int s = 0; s < spp; s++) {
						Spectrum L = f(scene, x, y, rng);
						if (isfinite(L)) {
							// Hacky: exclude NaNs in the rendering.
							radiance += L;
						}
					}
					img(x, y) = radiance / Real(spp);
				}
			}
			reporter.update(1);
		}, Vector2i(num_tiles_x, num_tiles_y));

		reporter.done();
	} else {

		ProgressReporter reporter(h*w);
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				if (!debug(x, y)) {
					continue;
				}
//			printf("x, y is %d, %d\n", x, y);
				Spectrum radiance = make_zero_spectrum();
				int spp = scene.options.samples_per_pixel;
				pcg32_state rng = init_pcg32(x*w + y);
				for (int s = 0; s < spp; s++) {
					if (debug(x, y)) {
						printf("\n\nsample id is %d\n", s);
						if (s == 2) {
							debug(x, y);
						}
					}
					Spectrum L = f(scene, x, y, rng);
					if (isfinite(L)) {
						// Hacky: exclude NaNs in the rendering.
						radiance += L;
					}
				}
				img(x, y) = radiance / Real(spp);
				if (debug(x, y)) {
					printf("Radiance for %d, %d is ", x, y);
					print(img(x, y));
				}
//				reporter.update(1);
			}

		}

		reporter.done();

	}


    return img;
}



Image3 photon_mapping_render(const Scene &scene) {
	int w = scene.camera.width, h = scene.camera.height;
	Image3 img(w, h);

	PhotonMapping photon_mapping(scene);
	photon_mapping.build_photon_map(scene);

	constexpr int tile_size = 16;
	int num_tiles_x = (w + tile_size - 1) / tile_size;
	int num_tiles_y = (h + tile_size - 1) / tile_size;

	if (0) {
		ProgressReporter reporter(h*w);
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				if (!debug(x, y)) {
					continue;
				}
				Spectrum radiance = make_zero_spectrum();
				int spp = scene.options.samples_per_pixel;
				pcg32_state rng = init_pcg32(x*w + y);
				for (int s = 0; s < spp; s++) {
					Spectrum L = photon_mapping.estimate_radiance(scene, x, y, rng);
					if (isfinite(L)) {
						// Hacky: exclude NaNs in the rendering.
						radiance += L;
					}
				}
				img(x, y) = radiance / Real(spp);
				if (debug(x, y)) {
					printf("Radiance for %d, %d is ", x, y);
					print(img(x, y));
				}
//				reporter.update(1);
			}
		}

		reporter.done();
	} else {
		ProgressReporter reporter(num_tiles_x * num_tiles_y);
		parallel_for([&](const Vector2i &tile) {
			// Use a different rng stream for each thread.
			pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
			int x0 = tile[0] * tile_size;
			int x1 = min(x0 + tile_size, w);
			int y0 = tile[1] * tile_size;
			int y1 = min(y0 + tile_size, h);
			for (int y = y0; y < y1; y++) {
				for (int x = x0; x < x1; x++) {
					Spectrum radiance = make_zero_spectrum();
					int spp = scene.options.samples_per_pixel;
					for (int s = 0; s < spp; s++) {
						radiance += photon_mapping.estimate_radiance(scene, x, y, rng);
					}
					img(x, y) = radiance / Real(spp);
				}
			}
			reporter.update(1);
		}, Vector2i(num_tiles_x, num_tiles_y));
		reporter.done();
	}
	return img;
}

Image3 render(const Scene &scene) {
    if (scene.options.integrator == Integrator::Depth ||
            scene.options.integrator == Integrator::ShadingNormal ||
            scene.options.integrator == Integrator::MeanCurvature ||
            scene.options.integrator == Integrator::RayDifferential ||
            scene.options.integrator == Integrator::MipmapLevel) {
        return aux_render(scene);
    } else if (scene.options.integrator == Integrator::Path) {
        return path_render(scene);
    } else if (scene.options.integrator == Integrator::VolPath) {
        return vol_path_render(scene);
	} else if (scene.options.integrator == Integrator::PhotonMapping) {
		return photon_mapping_render(scene);
	} else {
        assert(false);
        return Image3();
    }
}
