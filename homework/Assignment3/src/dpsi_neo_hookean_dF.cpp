#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
	dw.setZero();
	Eigen::Matrix3d dw33;
	//这样求伴随矩阵结果不对劲,太坑了，怎么这么坑，找了半天bug,结果这个库函数不行！
	//Eigen::Matrix3d Fadj = F.adjoint().transpose();
	Eigen::Matrix3d Fadj;
	Fadj << F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1), -F(1, 0)* F(2, 2) + F(1, 2) * F(2, 0), F(1, 0)* F(2, 1) - F(1, 1) * F(2, 0),
			-F(0, 1)* F(2, 2) + F(0, 2) * F(2, 1), F(0, 0)* F(2, 2) - F(0, 2) * F(2, 0), -F(0, 0)* F(2, 1) + F(0, 1) * F(2, 0),
			F(0, 1)* F(1, 2) - F(0, 2) * F(1, 1), -F(0, 0)* F(1, 2) + F(0, 2) * F(1, 0), F(0, 0)* F(1, 1) - F(0, 1) * F(1, 0);

    double J = F.determinant();
	double Jm23 = 1.0/pow(J, 2.0 / 3.0);
	double Jm53 = 1.0/pow(J, 5.0 / 3.0);

	Eigen::Matrix3d FT = F.transpose();
	Eigen::Matrix3d FT_F = FT * F;
	
	double tr = FT_F.trace();

	dw33 = C*Jm23*2.0*F + (-(2.0/3.0)*C*Jm53*tr+ D*2.0*(J-1.0))*Fadj;
	dw.segment(0,3) = dw33.col(0);
	dw.segment(3,3) = dw33.col(1);
	dw.segment(6,3) = dw33.col(2);

	return;
}