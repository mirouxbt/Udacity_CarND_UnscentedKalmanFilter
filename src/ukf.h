#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double _2PI = 2. * M_PI;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  
  ///* Radar noise measurement matrix
  MatrixXd R_radar_;

  ///* Radar noise measurement matrix
  MatrixXd R_lidar_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;
  
  ///* Number of sugma points
  int n_sigma_points_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage &meas_package);




private:
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage &meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage &meas_package);
  
  /**
   * GenerateSigmatPoints from the current state
   * @param Xsig_out must have the size (n_aug_, 2 * n_aug_ + 1)
   */
  void GenerateAugmentedSigmaPoints(MatrixXd &Xsig_out);
  
  /**
   * Predict the sigma points using the CTRV model
   * @param delta_t Time between k and k+1 in s
   * @param Xsig_in sigma points before prediction of size (n_aug_, 2 * n_aug_ + 1)
   * @param Xsig_out sigma point after prediction size (n_x_, 2 * n_aug_ + 1)
   */
  void PredictSigmaPoints(double delta_t, const MatrixXd &Xsig_in, MatrixXd &Xsig_out);
  
  /**
   * Keep the angle between -PI and PI
   * @param angle to keep to [-PI, PI]
   * @return angle normalized 
   */
  double AngleNormalization(double angle) {
    while (angle >  M_PI) angle -= _2PI;
    while (angle < -M_PI) angle += _2PI;
    return angle;
  }

  /**
   * Predict the sigma point using the Radar measurement
   * Will update the matrix vector and matrix
   * @param z_pred_out prediction measurement must have the size (3)
   * @param S_out      measurement covariance must have the size (3,3)
   * @param Zsig_out   sigma point in measurement space must hace the size (3, 2 * n_aug_ + 1)
   */  
  void PredictRadarMeasurement(VectorXd &z_pred_out, MatrixXd &S_out, MatrixXd &Zsig_out);

  /**
   * Predict the sigma point using the Lidar measurement
   * Will update the matrix vector and matrix
   * @param z_pred_out prediction measurement must have the size (3)
   * @param S_out      measurement covariance must have the size (3,3)
   * @param Zsig_out   sigma point in measurement space must hace the size (3, 2 * n_aug_ + 1)
   */  
  void PredictLidarMeasurement(VectorXd &z_pred_out, MatrixXd &S_out, MatrixXd &Zsig_out);

};

#endif /* UKF_H */
