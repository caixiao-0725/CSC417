//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    //这里写的2阶，4阶太麻烦了，放弃  
    Eigen::VectorXd f;
    Eigen::VectorXd k1;
    Eigen::VectorXd kdot1;
    
    force(f, q, qdot);
    k1 = f / mass;
    kdot1 = qdot;

    Eigen::VectorXd k2;
    Eigen::VectorXd kdot2;
    Eigen::VectorXd yt_add_tk1;
    Eigen::VectorXd yt_add_tkdot1;
    yt_add_tk1 = yt_add_tk1 = q + dt * k1;
    yt_add_tkdot1 = qdot + dt * kdot1;
    force(f, yt_add_tk1, yt_add_tkdot1);
    k2 = f / mass;
    kdot2 = yt_add_tkdot1;

    qdot += (dt / 2) * (k1 + k2);
    q += (dt / 2) * (kdot1 + kdot2);
    
}