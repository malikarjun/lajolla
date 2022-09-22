//
// Created by Mallikarjun Swamy on 2/5/22.
//

#include <fstream>
#include "util.h"

bool between(int val, int low, int high) {
	return low <= val && val <= high;
}

bool debug(int x, int y) {
	return x == 256 && y == 256;
//	return false;
	return true;


//	return between(x, 105, 272) && between(y, 107, 243);
}


void print(Vector3 vec, std::string str) {
	std::cout << str << " " << vec << std::endl;
}

void print(Vector2 vec, std::string str) {
    std::cout << std::endl << str << " " << vec << std::endl;
}

void print(Real val) {
	std::cout << val << std::endl;
}

void save_mat(Matrix4x4 mat, const std::string& filename) {
    // store data in row-major format
    size_t size = 4;
    std::vector<Real> mat_data;
    std::vector<size_t> shape = {size, size};
    for (int i = 0; i < shape[0]; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            mat_data.push_back(mat.data[i][j]);
        }
    }
    cnpy::npy_save(filename, mat_data.data(), shape);
}

void save_npy(Image4 img, const std::string& filename) {
    std::vector<Real> vec;
    size_t w = img.width, h = img.height, c = 4;


    for (int x = 0; x < img.width; ++x) {
        for (int y = 0; y < img.height; ++y) {
            for (int k = 0; k < c; ++k) {
                vec.push_back(img(x, y)[k]);
            }
        }
    }

    std::vector<size_t> shape = {h, w, c};
    cnpy::npy_save(filename, vec.data(), shape);
}

void save_npy(std::vector<Matrix4x4> mats, const std::string& filename) {
    std::vector<Real> vec;
    size_t x = mats.size(), y = 4, z = 4;
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            for (int k = 0; k < z; ++k) {
                vec.push_back(mats[i].data[j][k]);
            }
        }
    }
    std::vector<size_t> shape = {x, y, z};
    cnpy::npy_save(filename, vec.data(), shape);

}


void save_txt(std::vector<std::string> text, const std::string& filename) {
    std::ofstream ofs;

    ofs.open (filename);
    for(auto it : text) {
        ofs << it << "\n";
    }
    ofs.close();
}

std::string replace(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return str;
    str.replace(start_pos, from.length(), to);
    return str;
}