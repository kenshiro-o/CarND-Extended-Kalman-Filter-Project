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

  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.H_ = H_laser_;

  // Initialise process noise
  noiseAx = 9;
  noiseAy = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  std::string mesurementType = "";

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    ekf_.x_ = VectorXd(4);

    // Set the state covariance matrix to random values
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
      mesurementType = "RADAR";
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];

      float x = rho * cos(phi);
      float y = rho * sin(phi);
      ekf_.x_ << x, y, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      /**
      Initialize state.
      */
      mesurementType = "LASER";
    }

    cout << "State vector initialised: " << mesurementType << " " << ekf_.x_ << endl;

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
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
  mesurementType = measurement_pack.sensor_type_ == MeasurementPackage::LASER ? "LASER" : "RADAR";
  // cout << "** NOT FIRST MEASUREMENT: " << mesurementType << endl;

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  float dt_sq = dt * dt;
  // cout << "current_time=" << measurement_pack.timestamp_
  //      << " previou_time=" << previous_timestamp_
  //      << " time_diff=" << dt << endl;
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

  MatrixXd aAT = MatrixXd(2, 2);
  aAT << noiseAx, 0,
      0, noiseAy;

  // cout << "Calculating the process noise covariance matrix" << endl;
  // Setting the process noise covariance matrix Q
  ekf_.Q_ = G * aAT * G.transpose();

  ekf_.Predict();
  // cout << "x Vector prediction done: " << ekf_.x_ << " " << mesurementType << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  // cout << "Just before update: " << mesurementType << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar measurements contains three values
    VectorXd z = VectorXd(3);
    z << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1],
        measurement_pack.raw_measurements_[2];

    ekf_.R_ = R_radar_;

    // Adding the measurement noise coefficient to our measurements
    // VectorXd w = VectorXd(3);
    // w << ekf_.R_(0, 0), ekf_.R_(1, 1), ekf_.R_(2, 2);
    z += ekf_.R_.diagonal();

    ekf_.UpdateEKF(z);
  }
  else
  {
    // Our LIDAR can only detect objects in the 2D cartesian space
    VectorXd z = VectorXd(2);
    z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
    ekf_.R_ = R_laser_;
    z += ekf_.R_.diagonal();
    ekf_.Update(z);
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
