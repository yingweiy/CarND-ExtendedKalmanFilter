#include "FusionEKF.h"
#include "kalman_filter.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_ = KalmanFilter();

  //create a 4D state vector, we don't know yet the values of the x state
  VectorXd x_ = VectorXd(4);

  //state covariance matrix P
  MatrixXd P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;


  //measurement matrix
  H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;

   //the initial transition matrix F_
  MatrixXd F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  MatrixXd Q_ = MatrixXd(4, 4);

  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);

  tools = Tools();

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;



    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        double rho = measurement_pack.raw_measurements_(0);
        double theta = measurement_pack.raw_measurements_(1);
        double p_x = rho*cos(theta);
        double p_y = rho*sin(theta);
        ekf_.x_ << p_x, p_y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // modify F and Q
  // update delta_t
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //2. Set the process covariance matrix Q
  // set the acceleration noise components
  float noise_ax = 5;
  float noise_ay = 5;
  ekf_.Q_ = MatrixXd(4,4);
  float dt4 = pow(dt, 4.0) / 4.0;
  float dt3 = pow(dt, 3.0) / 2.0;
  float ax2 = noise_ax * noise_ax;
  float ay2 = noise_ay * noise_ay;
  ekf_.Q_ << dt4 * ax2, 0, dt3 * ax2, 0,
           0, dt4*ay2, 0, dt3 * ay2,
           dt3 * ax2, 0, dt*dt*ax2, 0,
           0, dt3*ay2, 0, dt*dt*ay2;
  cout << "pre prediction...";
  ekf_.Predict();
    cout << "post prediction...";
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.R_ = R_radar_;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
