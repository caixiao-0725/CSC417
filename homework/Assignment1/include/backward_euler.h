#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force, STIFFNESS &stiffness) {
	Eigen::MatrixXd k_matrix;
	stiffness(k_matrix, q, qdot);
	float k = - k_matrix(0, 0);
	float det_A = 1 + k * dt * dt / mass;
	Eigen::VectorXd q_t;
	Eigen::VectorXd qdot_t;
	q_t = q;
	qdot_t = qdot;
	qdot = (qdot_t - (k * dt / mass) * q_t) / det_A;
	q = (dt * qdot_t + q_t) / det_A;
}