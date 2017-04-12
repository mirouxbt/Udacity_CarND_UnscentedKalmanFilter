#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Num sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;
  
  // Init the weights
  weights_ = VectorXd(n_sigma_points_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_(1) = 0.5 / (lambda_ + n_aug_);
  
  for (int i=2; i<n_sigma_points_; i++) {
    weights_(i) = weights_(1);
  }

  //Measurement noise covariance matrix
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0                      , 0,
              0                  , std_radphi_*std_radphi_, 0,
              0                  , 0                      ,std_radrd_*std_radrd_
  ;

 
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_*std_laspx_ , 0 ,
              0                     , std_laspy_*std_laspy_
  ;
  
  // Use the process noise longitudinal acc to init the covariance matrix P
  P_ *= std_a_*std_a_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Init state as much as we can
      float ro     = meas_package.raw_measurements_(0);
      float phi    = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      
      if (fabs(ro) < 0.001) ro = 0.001;
      
      x_ <<  ro * cos(phi), ro * sin(phi),
             ro_dot , phi , 0.0
      ;  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      float px = meas_package.raw_measurements_(0);
      float py = meas_package.raw_measurements_(1);

      if (fabs(px) < 0.001) px = 0.001;
      if (fabs(py) < 0.001) py = 0.001;
      
      x_ <<  px , py ,
             0.0, 0.0, 0.0
      ;  
    }

    // init timestamp
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  if (   (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
      || (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) ) {  
  
    // Compute elapsed time since last measurement
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    Prediction(delta_t);
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      UpdateRadar(meas_package);
    } 
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Laser updates
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // Generate augmented sigma point to take into consideration the process noise
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_points_);
  GenerateAugmentedSigmaPoints(Xsig_aug);

 
  // Predict sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);
  PredictSigmaPoints(delta_t, Xsig_aug, Xsig_pred_);

  /**
   *  Predict state and covariance using the sigma point
   */
    
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = AngleNormalization(x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &meas_package) {
  // Dimension of measurement space
  int n_z = 2;

  // Predixt sigma points in measurment space
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S      = MatrixXd(n_z,n_z);
  MatrixXd Zsig   = MatrixXd(n_z, n_sigma_points_);
  PredictLidarMeasurement(z_pred, S, Zsig);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  
  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = AngleNormalization(z_diff(1));
  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = AngleNormalization(x_diff(3));
  
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  // S inverse
  MatrixXd Si = S.inverse();
  
  //Kalman gain K;
  MatrixXd K = Tc * Si;
  
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  
  //angle normalization
  z_diff(1) = AngleNormalization(z_diff(1));
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();
  
  /**
   * NIS 
   */
  NIS_laser_ = z_diff.transpose() * Si * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &meas_package) {
  // Dimension of measurement space
  int n_z = 3;

  // Predixt sigma points in measurment space
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S      = MatrixXd(n_z,n_z);
  MatrixXd Zsig   = MatrixXd(n_z, n_sigma_points_);
  PredictRadarMeasurement(z_pred, S, Zsig);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  
  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = AngleNormalization(z_diff(1));
  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = AngleNormalization(x_diff(3));
  
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // S inverse
  MatrixXd Si = S.inverse();
  
  //Kalman gain K;
  MatrixXd K = Tc * Si;
  
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  
  //angle normalization
  z_diff(1) = AngleNormalization(z_diff(1));
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();
  
  /**
   * NIS 
   */
  NIS_radar_ = z_diff.transpose() * Si * z_diff;
}

/**
 * GenerateSigmatPoints from the current state
 * Will update the matrix Xsig_out
 * Xsig_out must have the size (n_aug_, 2 * n_aug_ + 1)
 */
void UKF::GenerateAugmentedSigmaPoints(MatrixXd &Xsig_out) {
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  //create augmented mean state
  x_aug << x_ , 0.0, 0.0;
  
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_aug_-2, n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd Ps = P_aug.llt().matrixL();
  Ps *= sqrt(lambda_ + n_aug_);
  
  //create augmented sigma points
  Xsig_out.col(0) = x_aug;
  
  for (int i=1; i<=n_aug_; ++i) {
      Xsig_out.col(i)        = x_aug;
      Xsig_out.col(i)       += Ps.col(i-1);
  
      Xsig_out.col(i+n_aug_)  = x_aug;
      Xsig_out.col(i+n_aug_) -= Ps.col(i-1);
  }
  
}

/**
 * Predict the sigma point using the CTRV model
 * Will update the matrix Xsig_out
 * Xsig_out must have the size (n_x_, 2 * n_aug_ + 1)
 */
void UKF::PredictSigmaPoints(double delta_t, const MatrixXd &Xsig_in, MatrixXd &Xsig_out) {
  
  // Pre-compute some value to avoid multiple computation
  
  
  //predict sigma points
  for (int i = 0; i< n_sigma_points_; i++)
  {
    //extract values from sigma input
    double p_x      = Xsig_in(0,i);
    double p_y      = Xsig_in(1,i);
    double v        = Xsig_in(2,i);
    double yaw      = Xsig_in(3,i);
    double yawd     = Xsig_in(4,i);
    double nu_a     = Xsig_in(5,i);
    double nu_yawdd = Xsig_in(6,i);
    
    double cos_y = cos(yaw);
    double sin_y = sin(yaw);

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        p_x += v/yawd * ( sin(yaw + yawd*delta_t) - sin_y);
        p_y += v/yawd * ( cos_y   - cos(yaw+yawd*delta_t) );
    }
    else {
        p_x += v * delta_t * cos_y;
        p_y += v * delta_t * sin_y;
    }

    yaw += yawd * delta_t;

    //add noise
    double noise_a  = nu_a * delta_t;
    double noise_aa = noise_a * 0.5 * delta_t;
    double noise_y  = nu_yawdd * delta_t;
    double noise_yy = noise_y * 0.5 * delta_t;
    
    p_x  += noise_aa * cos_y;
    p_y  += noise_aa * sin_y;
    v    += noise_a;
    yaw  += noise_yy;
    yawd += noise_y;

    //write output
    Xsig_out(0,i) = p_x;
    Xsig_out(1,i) = p_y;
    Xsig_out(2,i) = v;
    Xsig_out(3,i) = yaw;
    Xsig_out(4,i) = yawd;
  }
  
}


/**
 * Predict the sigma point using the Radar measurement
 * Will update the matrix vector and matrix
 * z_pred_out must have the size (3)
 * S_out must have the size (3,3)
 */
void UKF::PredictRadarMeasurement(VectorXd &z_pred_out, MatrixXd &S_out, MatrixXd &Zsig_out) {
  
  //predict sigma points in measurements sapce
  for (int i = 0; i< n_sigma_points_; i++)
  {
    //extract values from sigma input
    double px      = Xsig_pred_(0,i);
    double py      = Xsig_pred_(1,i);
    double v       = Xsig_pred_(2,i);
    double yaw     = Xsig_pred_(3,i);
    
    double vcos = v * cos(yaw);
    double vsin = v * sin(yaw);
    double rho_norm = sqrt(px*px + py*py);

    Zsig_out(0,i) = rho_norm;                       // Rho
    Zsig_out(1,i) = atan2(py, px);                  // Phi
    Zsig_out(2,i) = (px*vcos + py*vsin) / rho_norm; // Rho_dot

  }
  
  //mean predicted measurement
  z_pred_out.fill(0.0);
  for (int i=0; i < n_sigma_points_; i++) {
      z_pred_out += weights_(i) * Zsig_out.col(i);
  }

  //measurement covariance matrix S
  S_out.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  
    //residual
    VectorXd z_diff = Zsig_out.col(i) - z_pred_out;

    //angle normalization
    z_diff(1) = AngleNormalization(z_diff(1));
    
    S_out += weights_(i) * z_diff * z_diff.transpose(); 
  }

  //add measurement noise covariance matrix
  S_out += R_radar_;

}

/**
 * Predict the sigma point using the Lidar measurement
 * Will update the matrix vector and matrix
 * z_pred_out must have the size (2)
 * S_out must have the size (2,2)
 */
void UKF::PredictLidarMeasurement(VectorXd &z_pred_out, MatrixXd &S_out, MatrixXd &Zsig_out) {
  
  //predict sigma points in measurements sapce
  for (int i = 0; i< n_sigma_points_; i++)
  {
    Zsig_out(0,i) = Xsig_pred_(0,i);
    Zsig_out(1,i) = Xsig_pred_(1,i);
  }
  
  //mean predicted measurement
  z_pred_out.fill(0.0);
  for (int i=0; i < n_sigma_points_; i++) {
      z_pred_out += weights_(i) * Zsig_out.col(i);
  }

  //measurement covariance matrix S
  S_out.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  
    //residual
    VectorXd z_diff = Zsig_out.col(i) - z_pred_out;

    //angle normalization
    z_diff(1) = AngleNormalization(z_diff(1));
    
    S_out += weights_(i) * z_diff * z_diff.transpose(); 
  }

  //add measurement noise covariance matrix
  S_out += R_lidar_;

}
