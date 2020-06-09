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


FpopInference fpop_analytic_inference_recycle(PiecewiseSquareLosses * cost_model_fwd,
        PiecewiseSquareLosses * cost_model_rev,
        double * data_vec, int data_count,
        double decay_rate,
        double penalty,
        int thj,
        int window_size,
        double sig) {
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
  check_selective_inference(&model, thj, window_size, data_count, data_vec, decay_rate, penalty, sig, verbose);
  double pval = calc_p_value(&model, thj, window_size, data_count, data_vec, decay_rate, sig, false, verbose);
  double approximation_error = 0.0;
  return FpopInference(pval, approximation_error, model, thj);
}