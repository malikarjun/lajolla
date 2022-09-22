import os
import argparse
from os.path import join, exists


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lajolla_exec", help="lajolla executable")
parser.add_argument("-s", "--scene", help="scene xml data")
parser.add_argument("-o", "--output_data_path", help="output data path")

args = parser.parse_args()

exec_file = args.lajolla_exec
xml_file = os.path.join(args.scene, "frame{}.xml")
out_dir = args.output_data_path

os.makedirs(out_dir, exist_ok=True)

buffer_names = ["frame{}.exr", "frame{}_depth.exr", "frame{}_normal.exr"]

# ground truth file
full_out_file = join(out_dir, "frame1_gt.exr")
cmd = "{exec_file} {xml_file} -o {out_file}".format(exec_file=exec_file, xml_file=xml_file.format("1_gt"),
                                                    out_file=full_out_file)

for frame in range(2):
    for name in buffer_names:
        full_out_file = join(out_dir, name.format(frame))
        cmd = "{exec_file} {xml_file} -o {out_file}".format(exec_file=exec_file, xml_file=xml_file.format(frame),
                                                            out_file=full_out_file)
        os.system(cmd)