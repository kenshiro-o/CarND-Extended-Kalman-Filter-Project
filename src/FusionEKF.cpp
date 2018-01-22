#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  // Declaration of Hj_, the Jacobian matrix use for RADAR measurements
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.H_ = H_laser_;

  // Initialisation of process noise covariance matrix
  // Both noise_ax and noise_ay are set to 9 for this exercise
  aAT_ = MatrixXd(2, 2);
  aAT_ << 9, 0,
      0, 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  std::string measurementType = "";

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    // first measurement
    ekf_.x_ = VectorXd(4);

    // Set the state covariance matrix to values from the Udacity lessons
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      measurementType = "RADAR";
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];

      float x = rho * cos(phi);
      float y = rho * sin(phi);
      ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      measurementType = "LASER";
    }

    cout << "State vector initialised: " << measurementType << " " << ekf_.x_ << endl;

    // Let's not forget to set the timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  measurementType = measurement_pack.sensor_type_ == MeasurementPackage::LASER ? "LASER" : "RADAR";

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  float dt_sq = dt * dt;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Reset the prediction coefficient matrix (only dt changes)
  ekf_.F_ << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  // Set the G matrix
  MatrixXd G = MatrixXd(4, 2);
  G << dt_sq / 2, 0,
      0, dt_sq / 2,
      dt, 0,
      0, dt;

  // Setting the process noise covariance matrix Q
  ekf_.Q_ = G * aAT_ * G.transpose();
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar measurements contains three values
    VectorXd z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1],
        measurement_pack.raw_measurements_[2];

    ekf_.R_ = R_radar_;
    // Adding the RADAR measurement noise coefficient to our measurements
    z += ekf_.R_.diagonal();

    ekf_.UpdateEKF(z);
  }
  else
  {
    // Our LIDAR can only detect objects in the 2D cartesian space
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];

    ekf_.R_ = R_laser_;
    // Adding the LASERmeasurement noise coefficient to our measurements
    z += ekf_.R_.diagonal();
    ekf_.Update(z);
  }
}
