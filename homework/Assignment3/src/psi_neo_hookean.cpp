#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    // psi = C * (J^(-2/3)(F.transpose() * F).trace() - 3) + D * (J - 1) * (J - 1);
    double J = F.determinant();
    psi = C * (pow(J, -2.0 / 3.0) * (F.transpose() * F).trace() - 3) + D * pow(J - 1, 2.0);
    return;
}