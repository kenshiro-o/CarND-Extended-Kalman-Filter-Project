#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  VectorXd z_pred = H_ * x_;
  // cout << "[LASER] z Vector calculated " << endl;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();

  MatrixXd S = H_ * P_ * Ht + R_;
  // cout << "[LASER] S Matrix calculated " << endl;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  // cout << "[LASER] PHt Matrix calculated " << endl;
  MatrixXd K = PHt * Si;
  // cout << "[LASER] K Matrix calculated " << endl;

  //new estimate
  x_ = x_ + (K * y);
  // cout << "[LASER] x Vector updated " << endl;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  // cout << "[LASER] P Matrix calculated " << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  // Define the individual variables from vector x for better readability
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // our h(xPrime) vector consist of
  // - range rho
  // - angle phi
  // -radial velocity rhoDot
  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rhoDot = (px * vx + py * vy) / rho;
  VectorXd h = VectorXd(3);
  h << rho, phi, rhoDot;

  // TODO adjust the phi angle so that it is between -pi and pi
  VectorXd y = z - h;
  phi = y(1);
  // cout << "phi=" << phi << endl;
  if (phi > M_PI)
  {
    // cout << "Phi angle greater than pi: " << phi;
    phi = phi - M_PI * 2.0;
    // cout << "Phi now in [-pi, pi]: " << phi;
  }
  else if (phi < -M_PI)
  {
    // cout << "Phi angle less than -pi: " << phi;
    phi += M_PI_2 * 2.0;
    // cout << "Phi now in [-pi, pi]: " << phi;
  }
  y(1) = phi;

  MatrixXd Hj = tools_.CalculateJacobian(x_);

  // cout << "Jacobian calculated " << endl;

  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_;
  // cout << "S Matrix calculated " << endl;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  // cout << "PHt Matrix calculated " << endl;
  MatrixXd K = PHt * Si;
  // cout << "K Matrix calculated " << endl;

  //new estimate
  x_ = x_ + (K * y);
  // cout << "x Vector updated" << endl;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
