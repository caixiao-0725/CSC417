#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {
    Eigen::Matrix36d J1;
    Eigen::Matrix66d J2;
    J1.setZero();
    J2.setZero();
    J1(0, 3) = 1;
    J1(1, 4) = 1;
    J1(2, 5) = 1;
    J1(0, 1) = -x(2);
    J1(0, 2) = x(1);
    J1(1, 0) = x(2);
    J1(1, 2) = -x(0);
    J1(2, 0) = -x(1);
    J1(2, 1) = x(0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            J2(i,j) = R(j,i);
        }
    }
    for (int i = 3; i < 6; i++) {
        J2(i, i) = 1;
    }
    J = R * J1 * J2;
   
}

