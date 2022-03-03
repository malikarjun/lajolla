//
// Created by Mallikarjun Swamy on 2/5/22.
//

#include "util.h"

bool between(int val, int low, int high) {
	return low <= val && val <= high;
}

bool debug(int x, int y) {
	return x == 360 && y == 5;
	return false;
	return true;


//	return between(x, 105, 272) && between(y, 107, 243);
}


void print(Vector3 vec, std::string str) {
	std::cout << str << " " << vec << std::endl;
}

void print(Real val) {
	std::cout << val << std::endl;
}