import os


exec_file = "/Users/mallikarjunswamy/imp/acads/courses/winter-2022/CSE_272/lajolla_public/cmake-build-debug/lajolla"
xml_file = "/Users/mallikarjunswamy/imp/acads/courses/winter-2022/CSE_272/lajolla_public/scenes/cbox/frame{}.xml"
out_file = "/Users/mallikarjunswamy/imp/acads/diffrender/svgf/pysvgf/data_fixed"

buffer_names = ["frame{}.exr", "frame{}_depth.exr", "frame{}_normal.exr"]

for frame in range(2):
    for name in buffer_names:
        full_out_file = os.path.join(out_file, name.format(frame))
        cmd = "{exec_file} {xml_file} -o {out_file}".format(exec_file=exec_file, xml_file=xml_file.format(frame),
                                                            out_file=full_out_file)
        os.system(cmd)