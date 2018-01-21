#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  // We have a vector rmse that will store the RMSE values for:
  // (rmse_pos_x, rmse_pos_y, rmse_v_x, rmse_v_y)
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  int estSize = estimations.size();
  int gtSize = ground_truth.size();

  // Size of estimations and ground truth vectors should be poisitive
  if (estSize == 0 || gtSize == 0)
  {
    cout << "Size of estimations or ground truth vector 0:" << std::endl;
    return rmse;
  }

  // They should also match...
  if (estSize != gtSize)
  {
    cout << "Size of estimations and ground truth vectors differ." << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for (int i = 0; i < estSize; ++i)
  {
    VectorXd est = estimations.at(i);
    VectorXd g = ground_truth.at(i);

    // cout << "Estimations: " << est << std::endl;
    // cout << "Ground Tuth: " << g << std::endl;

    VectorXd res = est - g;
    VectorXd res_pow = res.array().square();
    res_pow /= estSize;

    // cout << "Res Pow: " << res_pow
    //      << std::endl;
    rmse = rmse + res_pow;
  }

  // rmse /= estSize;
  rmse = rmse.array().sqrt();

  // cout << "RMSE = " << rmse << std::endl;

  return rmse;
  // return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
  MatrixXd Hj(3, 4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float sq_add = px * px + py * py;

  //check division by zero
  if (sq_add == 0.0)
  {
    cout << "Division by xero averted" << std::endl;
    return Hj;
  }

  // Pre compute common denominator terms
  float sqrt_denom = sqrt(sq_add);
  float pow_denom = pow(sq_add, 3.0 / 2.0);

  // TODO (vx * py - vy * px) could be stored as a variable 'vp' as used
  // either as 'vp' or '-vp' depending on context

  Hj << px / sqrt_denom, py / sqrt_denom, 0, 0,
      -py / sq_add, px / sq_add, 0, 0,
      (py * (vx * py - vy * px)) / pow_denom, (px * (vy * px - vx * py)) / pow_denom, px / sqrt_denom, py / sqrt_denom;

  return Hj;
}
