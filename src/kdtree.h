//
// Created by Mallikarjun Swamy on 3/5/22.
//

#ifndef LAJOLLA_KDTREE_H
#define LAJOLLA_KDTREE_H

#include <vector>
#include <algorithm>
#include <queue>
#include "vector.h"

class KDRecord
{
public:
	Vector3 position;
	Vector3 power;
	Vector3 dir_in;
	//TODO: Add Photon

	KDRecord() {}
	KDRecord(Vector3 position, Vector3 power, Vector3 dir_in) : position(position), power(power), dir_in(dir_in){}
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
		:left_node(nullptr),right_node(nullptr)
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
												 kdtree_closest_point(next_branch, heap, pivot, depth + 1,k),
												 pivot);
	if (!best.has_value())
		return best;

	if (heap.size() < k || powf((pivot[axis] - root->record.position[axis]), 2) < heap.top().dis * heap.top().dis)
	{
		best = closest_point(best,
							 kdtree_closest_point(opposite_branch, heap, pivot, depth + 1,k),
							 pivot);
	}

	return best;
}


std::vector<KDRecord> kdtree_knn(KDTree* root, Vector3 pivot, int k)
{

	std::priority_queue<KNNRecord, std::vector<KNNRecord>, KNNCompare> heap;
	kdtree_closest_point(root, heap, pivot, 0,k);
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
												 kdtree_closest_point(next_branch,pivot,depth+1),
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

Vector3 closest_test(const std::vector<Vector3>& points,Vector3 pivot)
{
	Vector3 best;
	Real best_dis = std::numeric_limits<Real>::max();
	for (int i = 0; i < points.size(); i++)
	{
		Real dis = distance(pivot, points[i]);
		if (dis  < best_dis)
		{
			best_dis = dis;
			best = points[i];
		}
	}
	return best;
}

std::vector<Vector3> knn_test(const std::vector<Vector3>& points, Vector3 pivot,int k)
{
	std::vector<std::pair<Real,Vector3>> dis_rec;
	for (int i = 0; i < points.size(); i++)
	{
		Real dis = distance(pivot, points[i]);
		dis_rec.push_back(std::make_pair(dis, points[i]));
	}
	std::sort(dis_rec.begin(), dis_rec.end(),
			  [](const std::pair<Real, Vector3>& lhs, const std::pair<Real, Vector3>& rhs)
			  {
				  return lhs.first < rhs.first;
			  });
	std::vector<Vector3> res;
	for (int i = 0; i < points.size(); i++)
	{
		res.push_back(dis_rec[i].second);
	}
	return res;
}

/*Example

	int correct = 0;
	int max_test = 1000;
	std::default_random_engine generator;
	generator.seed(42);
	std::uniform_real_distribution<Real> distribution(0.0, 1000.0);
	for (int test_count = 0; test_count < max_test; test_count++)
	{
		int n = 20000;
		std::vector<Vector3> points;
		for (int i = 0; i < n; i++)
		{
			points.push_back(Vector3(distribution(generator), distribution(generator), distribution(generator)));
		}
		int k = 100;
		KDTree* kdtree = new KDTree(std::vector<Vector3>(points), 0);
		Vector3 pivot(distribution(generator), distribution(generator), distribution(generator));
		std::vector<KDRecord> v1 = kdtree_knn(kdtree, pivot, k);
		delete kdtree;
		std::vector<Vector3>  v2 = knn_test(points, pivot,k);
		bool all_correct_flag = true;
		for (int j = 0; j < k; j++)
		{
			all_correct_flag &= (distance(v1[j].point, v2[j]) <= std::numeric_limits<Real>::epsilon());
		}
		if (all_correct_flag)
		{
			correct += 1;
		}

	}

	std::cout << correct << "/" << max_test << std::endl;

*/

#endif //LAJOLLA_KDTREE_H
