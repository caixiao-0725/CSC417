#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
	double alpha = 1.0;
	for (int i = 0; i < maxSteps; i++) {
		tmp_g.setZero();
		g(tmp_g, x0);
		if (tmp_g.norm() < 1e-8) return tmp_g.norm();

		tmp_H.setZero();
		H(tmp_H, x0);
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.compute(tmp_H);

		Eigen::VectorXd dx = solver.solve(tmp_g);

		alpha = 1.0;
		while((f(x0 - alpha*dx) > f(x0) - 1e-8*alpha*tmp_g.dot(dx))&&(alpha>1e-12))
			   alpha *= 0.5;

		x0 -= alpha*dx;
   }
   return tmp_g.norm();
}
