#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	P.resize(q_size - 3 * indices.size(), q_size);
	P.setZero();
	int count = 0;

	//Ҫ��P��ÿһ�ж���1��1Ҫ����Ų
	//P 1 0 0 0 0 0 0 0 0 
	//  0 1 0 0 0 0 0 0 0
	//  0 0 1 0 0 0 0 0 0
	//  0 0 0 0 0 0 1 0 0
	//  0 0 0 0 0 0 0 1 0
	//  0 0 0 0 0 0 0 0 1

	for (int i = 0; i < q_size - 3 * indices.size(); i++) {
		if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
		P.insert(i, i + 3.0 * count) = 1.0;
	}

}