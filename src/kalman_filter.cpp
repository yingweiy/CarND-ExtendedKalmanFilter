#include "kalman_filter.h"
#include <iostream>
#define PI 3.1415926

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;


// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
    double px = x_(0);
    double py = x_(1);
    double vx  = x_(2);
    double vy = x_(3);
    double rho = sqrt(px*px + py*py);
    double rho_dot;

    if (rho<0.001)
    {
        rho_dot = 0.0;
    }
    else
    {
        rho_dot = (px*vx+py*vy)/rho;
    }

  VectorXd z_pred = VectorXd(3);
  VectorXd z_norm = VectorXd(3);

  z_pred << rho, atan2(py, px), rho_dot;
  //z_norm << z(0), atan2(sin(phi), cos(phi)), z(2);
  //z_norm << z(0), phi, z(2);
  cout << z(1) << "," << z(2) << endl;
  VectorXd y = z - z_pred;
  double dy = y(1);
  if (dy>PI) {
        y(1) = dy - 2*PI;
    }
  if (dy<-PI) {
        y(1) = dy + 2*PI;
    }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
