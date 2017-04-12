#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0.0, 0.0, 0.0, 0.0;

  // Validation
  if(    estimations[0].size() != 4 
      || estimations.size() != ground_truth.size() ) {
    std::cerr <<  "CalculateRMSE : Invalid or mismatch size of parameters" << std::endl;
    return rmse;
  }
  
  // accumulate diff
  for(size_t i=0; i<estimations.size(); ++i) {
    VectorXd acc = estimations[i] - ground_truth[i];
    acc = acc.array() * acc.array();
    
    rmse += acc;
  }
  
  // Average it
  rmse /= estimations.size();
  
  // Square root it
  rmse = rmse.array().sqrt();
  
  return rmse;
}
