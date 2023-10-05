#include <d2psi_neo_hookean_dq2.h>

void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
	//我感到很疲倦，离家乡已经很远，害怕再也不能回到你身边```
	ddw.setZero();
	double J = F.determinant();
	Eigen::Matrix3d Fadj = F.adjoint();

	double Jm23 = pow(J, -2.0 / 3.0);
	double Jm53 = pow(J, -5.0 / 3.0);

	Eigen::Matrix3d FT = F.transpose();
	Eigen::Matrix3d FT_F = FT * F;

	double tr = FT_F.trace();
	//第一项
	Eigen::Matrix99d I99;
	I99.setIdentity();
	ddw = (2 * C * Jm23)*I99;
	//第二项
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Eigen::Matrix3d temp;
			temp = -(4.0 / 3.0) * C * Jm53 * Fadj(i, j) * F;
			//eigen不熟，只能这样一个一个装填
			ddw.block(0, 3 * j + i, 3, 1) += temp.col(0);
			ddw.block(3, 3 * j + i, 3, 1) += temp.col(1);
			ddw.block(6, 3 * j + i, 3, 1) += temp.col(2);
		}
	}
	//第三项
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Eigen::Matrix3d temp;
			temp = -(4.0 / 3.0) * C * Jm53 * Fadj * F(i,j);
			//eigen不熟，只能这样一个一个装填
			ddw.block(0, 3 * j + i, 3, 1) += temp.col(0);
			ddw.block(3, 3 * j + i, 3, 1) += temp.col(1);
			ddw.block(6, 3 * j + i, 3, 1) += temp.col(2);
		}
	}

	//第四项


}