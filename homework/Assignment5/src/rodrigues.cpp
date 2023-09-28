#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
	Eigen::Vector3d dir;
	double theta;
	theta = omega.norm();
	dir = omega.normalized();
	Eigen::Matrix3d skew;
	skew(0, 0) = 0;
	skew(0, 1) = -dir(2);
	skew(0, 2) = dir(1);
	skew(1, 0) = dir(2);
	skew(1, 1) = 0;
	skew(1, 2) = -dir(0);
	skew(2, 0) = -dir(1);
	skew(2, 1) = dir(0);
	skew(2, 2) = 0;
	R = Eigen::Matrix3d::Identity() + sin(theta) * skew + (1 - cos(theta)) * skew * skew;
}