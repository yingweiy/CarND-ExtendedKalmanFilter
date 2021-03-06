#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
    residue_ = VectorXd(4);
    residue_ << 0, 0, 0, 0;
    last_estimation_size_ = 0;
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    for(int i=last_estimation_size_; i < estimations.size(); ++i){
        VectorXd err = estimations[i] - ground_truth[i];
        err = err.array() * err.array();
        residue_ += err;
    }
    last_estimation_size_ = estimations.size();

    rmse = residue_ / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    MatrixXd Hj(3,4);
    Hj << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //check division by zero
    double pxy = px*px + py*py;
    double spxy = sqrt(pxy);
    double pxy15 = pxy * spxy;

    if (pxy<0.000001) {
        cout << "Warining: Dividing by zero in Jacobian calculation!" << endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj << px/spxy, py/spxy, 0, 0,
            -py/pxy, px/pxy, 0, 0,
            py*(vx*py - vy*px)/pxy15, px*(vy*px - vx*py)/pxy15, px/spxy, py/spxy;

    return Hj;

}
