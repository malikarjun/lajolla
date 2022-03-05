//
// Created by Mallikarjun Swamy on 3/4/22.
//

#ifndef LAJOLLA_PHOTON_MAP_H
#define LAJOLLA_PHOTON_MAP_H

#include "scene.h"
#include <queue>

struct Photon {
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
};


class PhotonMap {
private:
	std::vector<Photon> photons;
public:
	PhotonMap() = default;

	int get_num_photons() { return photons.size(); }
	Photon& get_ith_photon(int i) { return photons[i]; }

	void add_photon(Photon& photon) { photons.push_back(photon); }
	void set_photons(std::vector<Photon>& photons) {
		this->photons = photons;
	}

	void build() {
		// TODO: @ziyang build KD tree here
		// kdtree.setPoints(photons.data(), photons.size());
		// kdtree.buildTree();
	}


	/**
	 *
	 * @param p
	 * @param k
	 * @param max_dist update this variable with the maximum distance (squared value) between p and all the photons
	 * @return
	 */
	std::vector<int> queryKNearestPhotons(const Vector3 &p, int k,
										  Real& max_dist2) const {
		// TODO: @ziyang seach KD tree to return K nearest photons
		// return kdtree.searchKNearest(p, k, max_dist);
		std::vector<int> nn_photons;

		// max heap of size k
		std::priority_queue<DistPhoton, std::vector<DistPhoton>, Compare> pq;

		for (int i = 0; i < photons.size(); ++i) {
			Photon photon = photons[i];
			pq.push(DistPhoton{distance(photon.position, p), i});

			if (pq.size() > k) {
				pq.pop();
			}
		}
		max_dist2 = pq.top().dist2 * pq.top().dist2;

		while (!pq.empty()) {
			nn_photons.push_back(pq.top().photon_idx);
			pq.pop();
		}

		return nn_photons;
	}
};

#endif //LAJOLLA_PHOTON_MAP_H
