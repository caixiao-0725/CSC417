#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    double l = (q0 - q1).norm();
    Eigen::Matrix66d BtB;
    BtB << 1, 0, 0, -1, 0, 0,
        0, 1, 0, 0, -1, 0,
        0, 0, 1, 0, 0, -1,
        -1, 0, 0, 1, 0, 0,
        0, -1, 0, 0, 1, 0,
        0, 0, -1, 0, 0, 1;
    Eigen::Vector6d q;
    q << q0(0),
        q0(1),
        q0(2),
        q1(0),
        q1(1),
        q1(2);
    H = -stiffness * BtB + stiffness * l0 * BtB / l - stiffness * l0 * (BtB * q * (BtB * q).transpose()) / (l * l * l);
}