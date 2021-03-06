#include "fpop.h"
#include "fpop_inference.h"
#include "funPieceList.h"
#include "selective_inference.h"
#include "utils.h"
#include <vector>
#include <Rcpp.h>
#include <algorithm>    // std::min
#include <math.h>

using namespace Rcpp;

/*
 *
 * Utility functions
 *
 */

// utility function for fpop intervals
NumericMatrix convert_loss(PiecewiseSquareLoss *p, int s) {
  const int ncols = 8;
  const int n_pieces = p -> piece_list.size();
  NumericMatrix out(n_pieces, ncols);
  
  int row_i = 0;
  for (SquareLossPieceList::iterator it = p -> piece_list.begin(); it != p-> piece_list.end(); it++) {
    out(row_i, 0) = it -> Square;
    out(row_i, 1) = it -> Linear;
    out(row_i, 2) = it -> Constant;
    out(row_i, 3) = it -> min_mean;
    out(row_i, 4) = it -> max_mean;
    out(row_i, 5) = it -> prev_mean;
    out(row_i, 6) = it -> data_i + 1;
    out(row_i, 7) = s;
    row_i++;
  }
  return out;
}

/*
 *
 * Estimation functions
 *
 */

// [[Rcpp::export(name = ".fpop")]]
List fpop_interface2
        (NumericVector data,
         double gam,
         double penalty,
         double min_mean,
         double max_mean,
         NumericVector cost_mat_r,
         IntegerVector end_vec_r,
         NumericVector mean_vec_r,
         IntegerVector intervals_mat_r){


  double *data_ptr = data.begin();
  int data_count = data.size();
  double *cost_mat = cost_mat_r.begin();
  int *end_vec = end_vec_r.begin();
  double *mean_vec = mean_vec_r.begin();
  int *intervals_mat = intervals_mat_r.begin();

  PiecewiseSquareLosses out = fpop(data_ptr, data_count, gam, penalty, min_mean, max_mean);
  decode_fpop(out, data_count, cost_mat, end_vec, mean_vec, intervals_mat, gam);
  
  List pw_losses(out.size());
  for (int i = 0; i < out.size(); i++) {
    pw_losses[i] = convert_loss(&out[i], i + 1);
  }
  
  return pw_losses;

}


/*
 *
 * Inference functions
 *
 */

// [[Rcpp::export(name = ".fpop_inference")]]
List fpop_inference_interface_recycle
        (NumericVector data,
         double decay_rate,
         double penalty,
         int window_size,
         double sig,
         int return_dev = 0,
         bool return_ci = false,
         bool two_sided = false,
         double alpha = 0.05,
         double mu = 0,
         double lower_trunc = -1e6) {

  double *data_ptr = data.begin();
  int data_count = data.size();
  int verbose = 0;
  double pval;

// forward pass
  PiecewiseSquareLosses cost_model_fwd = fpop(data_ptr, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

// backward pass
  double *data_vec_rev = reverse_data(data_ptr, data_count);
  PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

  std::list<int> ll = extract_changepoints(cost_model_fwd, data_count); //extract changepoints
  std::list<int>::iterator it;

  const int ncols = 5;// changepoint + 1 (in R notation), pval, approximation error, LCB, UCB
  const int nrows =  ll.size() ;
  NumericMatrix out_mat(nrows, ncols);
  List phi_intervals(nrows);

  int row_i = 0;
  std::vector<int> ll_vec(ll.begin(), ll.end());
  for (int k = 0; k < ll_vec.size(); ++k){
    int it = ll_vec[k];
    try {
      if (it > 0) {
          double *v = construct_v(data_count, it, window_size, decay_rate);
          double vTy = construct_vTy(data_ptr, v, data_count, it, window_size);
          free(v);
          if (vTy >= lower_trunc) {
              FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, data_ptr,
                                                                  data_count, decay_rate, penalty, it,
                                                                  window_size, sig, return_ci, two_sided, alpha, mu,
                                                                  lower_trunc);
              out_mat(row_i, 0) = out.thj + 1;
              out_mat(row_i, 1) = out.pval;
              out_mat(row_i, 2) = out.approximation_error;
              if (return_ci) {
                  out_mat(row_i, 3) = out.confidence_interval[0];
                  out_mat(row_i, 4) = out.confidence_interval[1];
              }
              if (return_dev) {
                  phi_intervals[row_i] = convert_loss(&out.model, 0);
              } else {
                  phi_intervals[row_i] = nullptr;
              }
              row_i++;
          }
      }
    } catch (std::exception &ex) {
      forward_exception_to_r(ex);
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
    }
  }

  return List::create(out_mat, phi_intervals, row_i);

}