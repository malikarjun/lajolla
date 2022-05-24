//
// Created by Mallikarjun Swamy on 2/5/22.
//

#ifndef LAJOLLA_UTIL_H
#define LAJOLLA_UTIL_H

#include "spectrum.h"
#include "matrix.h"
#include "cnpy.h"
#include "image.h"

bool debug(int x, int y);


void print(Vector3 vec, std::string str="");
void print(Real val);

void save_mat(Matrix4x4 mat, const std::string& filename);

void save_npy(Image4 img, const std::string& filename);

void save_npy(std::vector<Matrix4x4> mats, const std::string& filename);

void save_txt(std::vector<std::string> text, const std::string& filename);

std::string replace(std::string str, const std::string& from, const std::string& to);

#endif //LAJOLLA_UTIL_H
