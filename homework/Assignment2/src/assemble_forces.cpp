#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    f.resize(q.size());
    f.setZero();
    for (int i = 0; i < E.rows(); i++) {
        Eigen::Vector3d q0 = q.segment<3>(3*E(i,0));
        Eigen::Vector3d q1 = q.segment<3>(3*E(i,1));
        double l0_ = l0(i);
        double l = (q0 - q1).norm();
        Eigen::Vector3d f0 = k * (l - l0_) * (1 / l) * (q0 - q1);
        f(E(i, 0) * 3) += f0(0);
        f(E(i, 0) * 3+1) += f0(1);
        f(E(i, 0) * 3+2) += f0(2);
        f(E(i, 1) * 3) -= f0(0);
        f(E(i, 1) * 3 + 1) -= f0(1);
        f(E(i, 1) * 3 + 2) -= f0(2);
    }    
    for (int i = 0; i < q.size()/3; i++) {
        f.segment<3>(3*i, 3) -= Eigen::Vector3d(0, 9.8f, 0);
    }
    f = -f;
};