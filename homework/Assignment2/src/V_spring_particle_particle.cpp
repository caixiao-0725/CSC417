#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    V = 0.0;
    Eigen::Vector3d delta = q0 - q1;
    double l = delta.norm();
    V = 0.5*stiffness*(l-l0)*(l-l0);
}