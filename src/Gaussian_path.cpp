//[[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string.h>
#include "functions.h"

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

//'@useDynLib MUSS
//'@import Rcpp
//'@import RcppEigen
// [[Rcpp::export]]
List Gaussian_path(const Eigen::MatrixXd Z, const Eigen::VectorXd y, const Eigen::VectorXd tauList,
                   const NumericVector spike_params, const float slab_param, const Eigen::VectorXd beta_init,
                   const bool sigma_update, const float sigma_init, const float theta_init,
                   const float a, const float b, const float omega, const float kappa,
                   const float tolerance, const int max_iter, const bool return_g){

  VectorXd beta = beta_init;

  int n = Z.rows();
  int p = Z.cols();
  int L = spike_params.length();

  float sigma, sigma_sqr, theta, diff, spike_param, g_k;
  MatrixXd hatX(n,p), hatXX(p,p), hatSigma(p,p), hatXX_D(p,p), beta_path(L,p);
  VectorXd beta_old;
  NumericVector p_k, pen_params, g_l, iter_nums(L), sigma_path(L), theta_path(L);
  List g_List = List(L);
  NumericVector thresholds(L);
  float c_2, wi, threshold;

  for (int l = 0; l < L; l++){
    sigma = sigma_init;
    theta = theta_init;
    if (return_g == true){
      g_l = 0;
    }
    spike_param = spike_params[l];
    for (int k = 0; k < max_iter; k++){

      // E step for X
      Rcpp::List postX = posteriorX(Z, y, beta, sigma, tauList);
      hatX = postX["hatX"];
      hatXX = postX["hatXX"];
      hatSigma = postX["hatSigma"];

      // E step for gamma
      Rcpp::NumericVector beta_NV(wrap(beta));
      p_k = theta / (theta + pow(slab_param/spike_param,0.5) * exp(-pow(beta_NV,2)/2*(1/spike_param-1/slab_param))*(1-theta));
      pen_params = (1-p_k)/spike_param + p_k/slab_param;


      ////  M-step ////

      //update beta
      beta_old = beta;
      hatXX_D = hatXX;
      for (int j = 0; j < p; j++){
        hatXX_D(j,j) += pow(sigma,2) * pen_params[j];
      }
      beta = hatXX_D.inverse() * hatX.transpose() * y;


      // update theta
      theta = (sum(p_k)+a-1) / (p+a+b-2);

      // update sigma
      if (sigma_update==true){
        sigma_sqr = 1/(n+omega+2)*(kappa*omega + (y.transpose() * y)(0) - 2* (y.transpose() * hatX * beta)(0)
                                     + (beta.transpose() * hatXX * beta)(0));
        sigma = pow(sigma_sqr,0.5);
      }

      //compute g
      if (return_g == true){
        g_k = comp_g(false, y, Z, tauList, hatSigma, hatX, hatXX, p_k,
                     spike_param, slab_param, omega, kappa, a, b, sigma, theta, beta);
        g_l.push_back(g_k);
      }


      //break or not
      diff = pow(((beta-beta_old).transpose() * (beta-beta_old))(0), 0.5);
      if ((diff <= tolerance) || (k==max_iter-1)){
        beta_path.row(l) = beta;
        sigma_path(l) = sigma;
        theta_path(l) = theta;
        if (return_g == true){
          g_l.erase(g_l.begin());
          g_List(l) = g_l;
        }
        iter_nums(l) = k+1;
        break;
      }
    }
    //compute threshold at each spike_param
    c_2 = slab_param/spike_param;
    wi = (1-theta)/theta;
    threshold = sigma * pow(2 * spike_param * log(wi* pow(c_2,0.5))*c_2/(c_2-1), 0.5);

    thresholds(l) = threshold;
  }

  // thresholded beta_output for Gaussian
  VectorXd beta_output = beta;
  for (int j = 0; j < p; j++){
    if (fabs(beta_output(j)) < thresholds[L-1]){
      beta_output(j) = 0;
    }
  }

  Rcpp::NumericMatrix beta_path_(wrap(beta_path));
  Rcpp::NumericVector beta_output_(wrap(beta_output));
  Rcpp::NumericVector beta_thresholds(wrap(thresholds));
  if (return_g == true){
    return List::create(
      _["beta_path"] = beta_path_,
      _["beta_output"] = beta_output_,
      _["beta_thresholds"] = beta_thresholds,
      _["sigma_path"] = sigma_path,
      _["theta_path"] = theta_path,
      _["g_List"] = g_List,
      _["iter_nums"] = iter_nums,
      _["spike_params"] = spike_params,
      _["slab_param"] = slab_param);
  }else{
    return List::create(
      _["beta_path"] = beta_path_,
      _["beta_output"] = beta_output_,
      _["beta_thresholds"] = beta_thresholds,
      _["sigma_path"] = sigma_path,
      _["theta_path"] = theta_path,
      _["iter_nums"] = iter_nums,
      _["spike_params"] = spike_params,
      _["slab_param"] = slab_param);
  }
}
