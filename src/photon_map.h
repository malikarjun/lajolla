//
// Created by Mallikarjun Swamy on 3/4/22.
//
#pragma once

#ifndef LAJOLLA_PHOTON_MAP_H
#define LAJOLLA_PHOTON_MAP_H

#include "scene.h"
#include <queue>
//#include "kdtree.h"

class KDRecord
{
public:
	Vector3 position;
	Vector3 power;
	Vector3 dir_in;
	//TODO: Add Photon

	KDRecord() {}
	KDRecord(Vector3 position, Vector3 power, Vector3 dir_in) : position(position), power(power), dir_in(dir_in) {}
};

class KNNRecord
{
public:
	KDRecord kdrecord;
	Real dis = std::numeric_limits<Real>::max();
};

struct KNNCompare
{
	bool operator()(const KNNRecord& lhs, const KNNRecord& rhs)
	{
		return lhs.dis < rhs.dis;
	}
};

class KDTree
{
public:
	KDTree* left_node;
	KDTree* right_node;
	KDRecord record;

	KDTree(std::vector<KDRecord>& records, int depth)
		:left_node(nullptr), right_node(nullptr)
	{
		if (records.size() == 0)
		{
			return;
		}

		int axis = depth % 3;
		std::sort(records.begin(), records.end(),
			[axis](const KDRecord& lhs, const KDRecord& rhs)
			{
				return lhs.position[axis] < rhs.position[axis];
			});
		record = records[records.size() / 2];
		std::size_t const half_size = records.size() / 2;
		std::vector<KDRecord> split_left(records.begin(), records.begin() + half_size);
		std::vector<KDRecord> split_right(records.begin() + half_size + 1, records.end());
		if (split_left.size() > 0)
		{
			left_node = new KDTree(split_left, depth + 1);
		}
		if (split_right.size() > 0)
		{
			right_node = new KDTree(split_right, depth + 1);
		}

	}
	~KDTree()
	{
		if (left_node != nullptr)
		{
			delete left_node;
		}

		if (right_node != nullptr)
		{
			delete right_node;
		}

	}

};

std::optional<KDRecord> closest_point(const std::optional<KDRecord>& node, const std::optional<KDRecord>& sub, const Vector3& pivot)
{
	if (!node.has_value() && !sub.has_value())
	{
		return std::nullopt;
	}
	else if (!node.has_value())
	{
		return sub;
	}
	else if (!sub.has_value())
	{
		return node;
	}

	Real dis1 = distance(node->position, pivot);
	Real dis2 = distance(sub->position, pivot);
	return dis1 < dis2 ? node : sub;
}


std::optional<KDRecord> kdtree_closest_point(KDTree* root, std::priority_queue<KNNRecord,
	std::vector<KNNRecord>, KNNCompare>& heap, Vector3 pivot, int depth, int k)
{
	if (root == nullptr)
	{
		return std::nullopt;
	}

	int axis = depth % 3;

	KDTree* next_branch = nullptr;
	KDTree* opposite_branch = nullptr;

	if (pivot[axis] < root->record.position[axis])
	{
		next_branch = root->left_node;
		opposite_branch = root->right_node;
	}
	else
	{
		next_branch = root->right_node;
		opposite_branch = root->left_node;
	}

	{
		KNNRecord knnr;
		knnr.dis = distance(root->record.position, pivot);
		knnr.kdrecord = root->record;
		if (heap.size() < k)
		{
			heap.push(knnr);
		}
		else if (knnr.dis < heap.top().dis)
		{
			heap.pop();
			heap.push(knnr);
		}

	}

	std::optional<KDRecord> best = closest_point(std::optional<KDRecord>(root->record),
		kdtree_closest_point(next_branch, heap, pivot, depth + 1, k),
		pivot);
	if (!best.has_value())
		return best;

	if (heap.size() < k || powf((pivot[axis] - root->record.position[axis]), 2) < heap.top().dis * heap.top().dis)
	{
		best = closest_point(best,
			kdtree_closest_point(opposite_branch, heap, pivot, depth + 1, k),
			pivot);
	}

	return best;
}


std::vector<KDRecord> kdtree_knn(KDTree* root, Vector3 pivot, int k)
{

	std::priority_queue<KNNRecord, std::vector<KNNRecord>, KNNCompare> heap;
	kdtree_closest_point(root, heap, pivot, 0, k);
	std::vector<KDRecord> result;
	for (int i = 0; i < k; i++)
	{
		result.push_back(heap.top().kdrecord);
		heap.pop();
	}
	std::reverse(result.begin(), result.end());
	return result;
}


std::optional<KDRecord> kdtree_closest_point(KDTree* root, Vector3 pivot, int depth)
{
	if (root == nullptr)
	{
		return std::nullopt;
	}

	int axis = depth % 3;

	KDTree* next_branch = nullptr;
	KDTree* opposite_branch = nullptr;

	if (pivot[axis] < root->record.position[axis])
	{
		next_branch = root->left_node;
		opposite_branch = root->right_node;
	}
	else
	{
		next_branch = root->right_node;
		opposite_branch = root->left_node;
	}

	std::optional<KDRecord> best = closest_point(std::optional<KDRecord>(root->record),
		kdtree_closest_point(next_branch, pivot, depth + 1),
		pivot);
	if (!best.has_value())
		return best;

	if (distance_squared(pivot, best->position) > powf((pivot[axis] - root->record.position[axis]), 2))
	{
		best = closest_point(best,
			kdtree_closest_point(opposite_branch, pivot, depth + 1),
			pivot);
	}

	return best;
}

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
