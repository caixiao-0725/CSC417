#include <d2psi_neo_hookean_dq2.h>
void cross_matrix(Eigen::Matrix3d &res, Eigen::Vector3d v) {
	res<< 0, -v(2), v(1),
		v(2), 0, -v(0),
		-v(1), v(0), 0;
}

void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
	//我感到很疲倦，离家乡已经很远，害怕再也不能回到你身边```
	ddw.setZero();
	double J = F.determinant();
    Eigen::Matrix3d Fadj;
    Fadj << F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1), -F(1, 0) * F(2, 2) + F(1, 2) * F(2, 0), F(1, 0)* F(2, 1) - F(1, 1) * F(2, 0),
        -F(0, 1) * F(2, 2) + F(0, 2) * F(2, 1), F(0, 0)* F(2, 2) - F(0, 2) * F(2, 0), -F(0, 0) * F(2, 1) + F(0, 1) * F(2, 0),
        F(0, 1)* F(1, 2) - F(0, 2) * F(1, 1), -F(0, 0) * F(1, 2) + F(0, 2) * F(1, 0), F(0, 0)* F(1, 1) - F(0, 1) * F(1, 0);

	double Jm23 = pow(J, -2.0 / 3.0);
	double Jm53 = pow(J, -5.0 / 3.0);
	double Jm83 = pow(J, -8.0 / 3.0);

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
	Eigen::Vector9d dJ_dF;
	dJ_dF.setZero();
	dJ_dF.segment(0, 3) = Fadj.col(0);
	dJ_dF.segment(3, 3) = Fadj.col(1);
	dJ_dF.segment(6, 3) = Fadj.col(2);
	Eigen::Matrix99d res;
	res = dJ_dF * dJ_dF.transpose();
	ddw += (2 * D +((10.0/9.0)*C*Jm83*tr))* res;

	//第五项 根据2020的siggraph course 41页  https://graphics.pixar.com/library/DynamicDeformablesSiggraph2020/paper.pdf
	double alpha = 2 * D * (J - 1) - 2.0 / 3.0 * C * Jm53 * tr;
	Eigen::Vector3d f0, f1, f2;
	f0 = F.col(0);
	f1 = F.col(1);
	f2 = F.col(2);
	Eigen::Matrix3d cross_f0, cross_f1, cross_f2;
	cross_matrix(cross_f0, f0);
	cross_matrix(cross_f1, f1);
	cross_matrix(cross_f2, f2);
	res.setZero();
	res.block(0, 3, 3, 3) = -cross_f2;
	res.block(0, 6, 3, 3) = cross_f1 ;
	res.block(3, 0, 3, 3) = cross_f2 ;
	res.block(3, 6, 3, 3) = -cross_f0;
	res.block(6, 0, 3, 3) = -cross_f1;
	res.block(6, 3, 3, 3) = cross_f0 ;

	ddw += alpha* res;

    
	return;
}