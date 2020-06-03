/* -*- compile-command: "R CMD INSTALL .." -*- */

#include <vector>
#include <stdio.h>
#include "funPieceList.h"
#include <math.h>
#include <stdexcept>
#include <stdlib.h>


void check_min_operation(PiecewiseSquareLoss *cost, PiecewiseSquareLoss *min_prev_cost, PiecewiseSquareLoss *cost_prev, double penalty, int data_i) {
  int status = cost -> check_min_of(min_prev_cost, cost_prev);
  try {
    if(status){
      printf("Lambda = %.20e\n", penalty);
      printf("Error at data_i = %d status = %d\n", data_i, status);
      cost -> set_to_min_env_of(min_prev_cost, cost_prev, true);
      printf("=min_prev_cost\n");
      min_prev_cost -> print();
      printf("=cost_prev + %f\n", penalty);
      cost_prev -> print();
      printf("=new cost model\n");
      cost -> print();
      throw std::runtime_error("min(f, g) error. Please report!");
    }
  } catch(int e) {
    printf("An exception occured %d \n", e);
    throw std::runtime_error("min(f, g) error. Please report!");
  }
}

PiecewiseSquareLoss * fpop_update_i(PiecewiseSquareLoss *cost, PiecewiseSquareLoss *cost_prev, double data, double decay_rate, double penalty, int data_i, bool fwd, int verbose) {

  PiecewiseSquareLoss min_prev_cost;
  PiecewiseSquareLoss scaled_cost;

  min_prev_cost.set_to_unconstrained_min_of(cost_prev, verbose);
  min_prev_cost.set_prev_seg_end(data_i - 1);
  min_prev_cost.add(0, 0, penalty);

  if (fwd) {
    scaled_cost.set_to_scaled_fwd(cost_prev, decay_rate);
  } else {
    scaled_cost.set_to_scaled_rev(cost_prev, decay_rate);
  }
  cost -> set_to_min_env_of(&min_prev_cost, &scaled_cost, verbose);
  check_min_operation(cost, &min_prev_cost, &scaled_cost, penalty, data_i);
  cost -> add(0.5, - data, data * data / 2);
  return(cost);
}


PiecewiseSquareLosses fpop
        (double *data_vec, int data_count, double decay_rate, double penalty, double min_mean, double max_mean) {

  PiecewiseSquareLoss *cost, *cost_prev;
  PiecewiseSquareLoss min_prev_cost;
  PiecewiseSquareLosses cost_model_mat(data_count);

  int verbose=0;
  for(int data_i=0; data_i<data_count; data_i++){
    cost = &cost_model_mat[data_i];
    if(data_i==0){
      cost -> piece_list.emplace_back(0.5, - data_vec[0], data_vec[0] * data_vec[0] / 2, min_mean, max_mean, -1, false);
      cost_prev = cost;
    }else {
      cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, data_vec[data_i], decay_rate, penalty, data_i, 1, verbose);
    }
  }
  return cost_model_mat;
}

PiecewiseSquareLosses fpop_custom
        (double *data_vec, int data_count,
         PiecewiseSquareLoss *start_cost,
         double decay_rate,
         double penalty,
         bool fwd,
         int verbose
        ){

  PiecewiseSquareLoss *cost_prev;
  PiecewiseSquareLosses cost_model_mat(data_count);

  cost_prev = start_cost;

  // build cost functions for each data point, starting with start_cost cost function
  for(int data_i=0; data_i < data_count; data_i++){
    cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, data_vec[data_i], decay_rate, penalty, data_i, fwd, verbose);
  }
  return cost_model_mat;
}

void decode_fpop(PiecewiseSquareLosses cost_model_mat,
                 int data_count,
                 double *cost_mat,
                 int *end_vec,
                 double *mean_vec,
                 int *intervals_mat,
                 double decay_rate) {

  PiecewiseSquareLoss *cost;
  double best_cost, best_mean, prev_mean;
  int prev_seg_end=data_count;

  for(int i=0; i< data_count; i++){
    cost = &cost_model_mat[i];
    intervals_mat[i] = cost->piece_list.size();
    cost->Minimize
            (&best_cost, &best_mean,
             &prev_seg_end, &prev_mean);

    cost_mat[i] = best_cost;
  }

  // first step
  cost = &cost_model_mat[data_count - 1];
  cost->Minimize
          (&best_cost, &best_mean,
           &prev_seg_end, &prev_mean);


  int prev_seg_old = data_count - 1;
  int out_i=0;

  // loop over all prev. changepoints
  while(prev_seg_old >= 0){
    if (prev_seg_old < data_count - 1) {
      cost = &cost_model_mat[prev_seg_end];
      cost->Minimize
              (&best_cost, &best_mean,
               &prev_seg_end, &prev_mean);
    }
    double last_mean = best_mean * decay_rate;
    for (int t = prev_seg_old; t > prev_seg_end; t--){
      mean_vec[out_i] = last_mean / decay_rate;
      last_mean = last_mean / decay_rate;
      end_vec[out_i] = prev_seg_end;
      out_i++;
    }
    prev_seg_old = prev_seg_end;
  }
}

std::list<int> extract_changepoints(PiecewiseSquareLosses cost_model_mat, int data_count) {
  int prev_seg_old = data_count - 1;
  std::list<int> changepoints;
  // loop over all prev. changepoints
  while(prev_seg_old >= 0){
    prev_seg_old = cost_model_mat[prev_seg_old].getMinIntervalInd();
    changepoints.emplace_front(prev_seg_old);
  }
  return changepoints;
}


PiecewiseBiSquareLosses fpop_2d_custom_start(double *data_vec, int data_count,
                                             PiecewiseSquareLoss *start_cost,
                                             double penalty,
                                             double decay_rate,
                                             bool forward,
                                             double nuTy,
                                             double norm_constant,
                                             int verbose) {

  PiecewiseBiSquareLosses collection;
  PiecewiseBiSquareLoss start;
  PiecewiseBiSquareLoss min_over_mu;

  start.set_to_pw_u(start_cost);
  collection.piece_list.clear();
  collection.piece_list.emplace_back(start);

  // print out collections
  //TODO: remove after debug
  verbose = 0;

  if(verbose) {
    collection.print();
  }


  // build cost functions for each data point, starting with start_cost cost function
  for(int data_i=0; data_i < data_count; data_i++){
    min_over_mu = collection.min_u();

    if (verbose) {
      printf("--------\n");
      printf("Collection at data_i=%d\n", data_i);
      printf("Min over mu of all elements in collection at t-1\n");
      min_over_mu.print();
    }

    min_over_mu.add(0, 0, 0, 0, 0, penalty);
    collection.collect(&min_over_mu);

    // TODO: check
    // basically the effective window size should be equal to the length of this subset array
    int effective_window_size = data_count;

    if (forward) {
      collection.rescale(forward, decay_rate);

      if (verbose) {
        printf("after rescaling\n");
        collection.print();
      }

      double decay_rate2 = decay_rate * decay_rate;
      double c1 = decay_rate*(decay_rate2-1)/(decay_rate2-pow(decay_rate,(-2*(effective_window_size-1)))) / norm_constant;

//      printf("norm_constant %f \n", norm_constant);
//      printf("nuTy %f \n", nuTy);
//      printf("effective_window_size %i \n", effective_window_size);

//      printf("uphi coeff %f \n", c1*pow(decay_rate, effective_window_size-data_i-1));
//      printf("c1 coeff %f \n", c1);
//      printf("effective_window_size %i \n", effective_window_size);
//      printf("phi^2 coeff %f\n",0.5*pow(c1,2)*pow(decay_rate, 2 * (effective_window_size-data_i-1)));

      collection.add(0.5,
                     (-c1*pow(decay_rate, -1*(effective_window_size-data_i-1))*nuTy) - data_vec[data_i],
                     c1*pow(decay_rate, -1*(effective_window_size-data_i-1)),
                     0.5*pow(c1,2)*pow(decay_rate, -2 * (effective_window_size-data_i-1)),
                     -pow(c1,2)*pow(decay_rate, -2 * (effective_window_size-data_i-1))*nuTy -
                     c1*pow(decay_rate, -1*(effective_window_size-data_i-1))*data_vec[data_i],
                     0.5*pow(data_vec[data_i],2) + 0.5*pow(c1,2)*
                     pow(decay_rate, -2 * (effective_window_size-data_i-1))*pow(nuTy,2) +
                     c1*pow(decay_rate, -1*(effective_window_size-data_i-1))*nuTy*data_vec[data_i]
      );
    } else {
      collection.rescale(forward, decay_rate);
      double decay_rate2 = decay_rate * decay_rate;
      double c1 = (decay_rate2 - 1)/( pow(decay_rate, 2 * (effective_window_size)) -1) / norm_constant;

      collection.add(0.5,
                     (c1*pow(decay_rate, effective_window_size - data_i - 1)*nuTy) - data_vec[data_i],
                     -c1*pow(decay_rate, effective_window_size - data_i - 1),
                     0.5*pow(c1,2)*pow(decay_rate, 2 * (effective_window_size - data_i - 1)),
                     -pow(c1,2)*pow(decay_rate, 2 * (effective_window_size - data_i -1))*nuTy +
                     c1*pow(decay_rate, effective_window_size - data_i - 1)*data_vec[data_i],
                     0.5*pow(data_vec[data_i],2) + 0.5*pow(c1,2)*
                     pow(decay_rate, 2*(effective_window_size - data_i - 1))*
                     pow(nuTy,2) -
                     c1*pow(decay_rate, effective_window_size - data_i - 1)*nuTy*data_vec[data_i]);
    }

    if (verbose) {
      collection.print();
    }

  }
  return collection;
}