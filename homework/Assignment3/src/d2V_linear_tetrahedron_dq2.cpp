#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
      //Code to compute non-integrated hessian matrix goes here
       //和前面的一阶导差不多
       Eigen::Matrix34d x;
       x << q.segment<3>(element(0) * 3),
           q.segment<3>(element(1) * 3),
           q.segment<3>(element(2) * 3),
           q.segment<3>(element(3) * 3);
       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);

       Eigen::Matrix99d d2psi;
       d2psi_neo_hookean_dF2(d2psi, x * dphi, C, D);

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

       // For row major d2psi
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

       dV = B.transpose() * d2psi * B;
     
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    //或许这里可以考虑把所有3*3的小矩阵都正定化？速度可能比把一个12*12的大的矩阵正定化要快，而且效果会更好
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
