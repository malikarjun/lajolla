//
// Created by Mallikarjun Swamy on 3/4/22.
//

#ifndef LAJOLLA_PHOTON_MAP_H
#define LAJOLLA_PHOTON_MAP_H

#include "scene.h"
#include <queue>
#include "kdtree.h"


/*struct Photon {
	Vector3 power;
	Vector3 position;
	Vector3 dir_in;

	Photon() {}
	Photon(const Vector3& power, const Vector3& position, const Vector3& dir_in)
		: power(power), position(position), dir_in(dir_in) {}
};

// used for querying k nearest neighbours
struct DistPhoton {
	Real dist2;
	int photon_idx;

	DistPhoton() {}
	DistPhoton(const Real &dist2, const int &photon_idx): dist2(dist2), photon_idx(photon_idx) {}
};

class Compare
{
public:
	bool operator() (DistPhoton dp1, DistPhoton dp2) {
		return dp1.dist2 < dp2.dist2;
	}
};*/

using Photon = KDRecord;

class PhotonMap {
private:
	std::vector<Photon> photons;
	KDTree* kdtree;
public:
	PhotonMap(std::vector<Photon> photons) : photons(photons) {
		kdtree = new KDTree(photons, 0);
	}

	int get_num_photons() { return photons.size(); }
	Photon& get_ith_photon(int i) { return photons[i]; }

	/**
	 *
	 * @param p
	 * @param k
	 * @param max_dist update this variable with the maximum distance (squared value) between p and all the photons
	 * @return
	 */
	std::vector<Photon> queryKNearestPhotons(const Vector3 &p, int k,
										  Real& max_dist2) const {
		std::vector<Photon> knn_photons = kdtree_knn( kdtree, p, k);

		/*		Real min_dist = distance(p, _knn_photons[0].position);
		std::vector<Photon> knn_photons;
		for(Photon photon : _knn_photons) {
			// TODO: remove the hardcoded threshold and use something scene dependent.
			if (distance(p, photon.position) > 100 * min_dist) {
//				printf("breaking early!!");
//				break;
			}
			knn_photons.push_back(photon);
		}*/

		max_dist2 = distance_squared(p, knn_photons[knn_photons.size()-1].position);
		return knn_photons;
	}
};

#endif //LAJOLLA_PHOTON_MAP_H
