#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
	dw.setZero();
	Eigen::Matrix3d dw33;
	Eigen::Matrix3d Fadj = F.adjoint();

    double J = F.determinant();
	double Jm23 = pow(J, -2.0 / 3.0);
	double Jm53 = pow(J, -5.0 / 3.0);


	Eigen::Matrix3d FT = F.transpose();
	Eigen::Matrix3d FT_F = FT * F;
	
	double tr = FT_F.trace();

	dw33 = C*Jm23*2*F - (2/3)*C*Jm53*tr*Fadj + D*(2*(J-1))*Fadj;
	dw.segment(0,3) = dw33.col(0);
	dw.segment(3,3) = dw33.col(1);
	dw.segment(6,3) = dw33.col(2);
	
	return;
}