#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
	VectorXd y = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	//MatrixXd PHt = P_ * Ht;
	//MatrixXd S = H_ * PHt + R_;
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K =  P_ * Ht * Si;

	//new estimate
	x_ = x_ + (K * y);
	//long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(4,4);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

	float range = sqrt(px*px + py*py);

	//check division by zero
	//if (fabs(px) < 0.0001 || fabs(range) < 0.0001) {
	  //cout << "UpdateEKF () - Error - Division by Zero" << endl;
	  //return;
	//}

	float angle = atan2(py,px);
	float range_rate = (px*vx + py*vy) / range;

	VectorXd z_pred = VectorXd(3);
	z_pred << range, angle, range_rate;

	VectorXd y = z - z_pred;
	if (y(1) > 3.14)
	{
	    y(1) = y(1) - (float)(2.0 * 3.14);
	}
	if (y(1) < -3.14)
	{
	    y(1) = y(1) + (float)(2.0 * 3.14);
	}
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	//long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(4,4);
	P_ = (I - K * H_) * P_;
}
