#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    double l = (q0 - q1).norm();
    Eigen::Vector3d f0 = stiffness * (l - l0) * (1 / l) * (q0 - q1);
    f(0) = f0(0);
    f(1) = f0(1);
    f(2) = f0(2);
    f(3) = -f0(0);
    f(4) = -f0(1);
    f(5) = -f0(2);
}