#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       Eigen::Matrix34d x;
       x << q.segment<3>(element(0) * 3),
           q.segment<3>(element(1) * 3),
           q.segment<3>(element(2) * 3),
           q.segment<3>(element(3) * 3);
       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);
       Eigen::Vector9d dpsi;
       dpsi_neo_hookean_dF(dpsi, x * dphi, C, D);
       //printf("%f\n", dpsi(1));
       //这个B矩阵的含义，本质上相当于先算dV的下面9个值，再对他们求和乘-1 ,算出上面的3个值

       Eigen::MatrixXd B(9, 12);
       B.setZero();
       // Since dpsi is computed as column major
       //      dx/dX   D00        D10        D20        D30       
       //      dy/dX       D00        D10        D20        D30   
       //      dz/dX          D00        D10        D20        D30
       //      dx/dY   D01        D11        D21        D31       
       // B =  dy/dY =     D01        D11        D21        D31   
       //      dz/dY          D01        D11        D21        D31
       //      dx/dZ   D02        D12        D22        D32       
       //      dy/dZ       D02        D12        D22        D32   
       //      dz/dZ          D02        D12        D22        D32
       for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 4; j++) {
               B.block<3, 3>(i * 3, j * 3) = Eigen::Matrix3d::Identity() * dphi(j, i);
           }
       }

       // For row major dpsi
       //      dx/dX   D00        D10        D20        D30       
       //      dx/dY   D01        D11        D21        D31       
       //      dx/dZ   D02        D12        D22        D32       
       //      dy/dX       D00        D10        D20        D30   
       // B =  dy/dY =     D01        D11        D21        D31   
       //      dy/dZ       D02        D12        D22        D32   
       //      dz/dX          D00        D10        D20        D30
       //      dz/dY          D01        D11        D21        D31
       //      dz/dZ          D02        D12        D22        D32
       // for(int i = 0; i < 3; i++) {
       //     for(int j = 0; j < 3; j++) {
       //         B.block<3,1>(i, 3*j+i) = dphi.row(j).transpose();
       //     }
       // }
       dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}