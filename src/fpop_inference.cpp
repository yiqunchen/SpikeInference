#include "fpop.h"
#include "funPieceList.h"
#include "selective_inference.h"
#include "utils.h"
#include "fpop_inference.h"
#include <math.h>
#include <cstdio>
#include <iostream>

FpopInference::FpopInference(double p, double a, PiecewiseSquareLoss m, int t) {
  pval = p;
  approximation_error = a;
  model = m;
  thj = t;
}

FpopInference::FpopInference(double p, double a, PiecewiseSquareLoss m, int t, std::vector<double> ci) {
    pval = p;
    approximation_error = a;
    model = m;
    thj = t;
    confidence_interval = ci;
}


FpopInference fpop_analytic_inference_recycle(PiecewiseSquareLosses * cost_model_fwd,
        PiecewiseSquareLosses * cost_model_rev,
        double * data_vec, int data_count,
        double decay_rate,
        double penalty,
        int thj,
        int window_size,
        double sig,
        bool return_ci = false,
        bool two_sided = false,
        double alpha = 0.05,
        double mu = 0,
        double lower_trunc = -INFINITY){ // return CI defaults to false
  int verbose = 0;
  double *data_vec_rev = reverse_data(data_vec, data_count);
  PiecewiseSquareLoss model = thj_in_model(
          cost_model_fwd,
          cost_model_rev,
          thj, // changepoint of interest
          window_size, // size of window around thj
          data_count, // number of data points
          data_vec, // original data
          data_vec_rev, // flipped data
          decay_rate, // AR1 decay parameter
          penalty, // tuning parameter to penalize the number of spikes
          verbose);


  //check_selective_inference(&model, thj, window_size, data_count, data_vec, decay_rate, penalty, sig, verbose);

  // impose possible truncation by vTy:

  //free(it);
  double pval = calc_p_value(&model, thj, window_size, data_count, data_vec, decay_rate, sig, two_sided, mu, lower_trunc);
  free(data_vec_rev); // free the data
  double approximation_error = 0.0;
  if (return_ci){
      std::vector<double> thj_CI = compute_CI(&model, thj,  window_size, data_count, data_vec, decay_rate, sig, alpha,
                                              lower_trunc);
      //free(data_vec);
      if((thj_CI[1]<=thj_CI[0])|(fabs(thj_CI[1]-thj_CI[0])<1e-6)){
          printf("Numerical problems encountered when computing confidence intervals, do NOT trust the results! \n");
      }
      return FpopInference(pval, approximation_error, model, thj, thj_CI);
  }else{
      return FpopInference(pval, approximation_error, model, thj);
  }

}