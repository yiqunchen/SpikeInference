#include "selective_inference.h"
#include "fpop.h"
#include "funPieceList.h"
#include <vector>
#include "utils.h"
#include <algorithm>    // std::min
#include <math.h>
#include <stdexcept>

/*
 * Contrast vector v
 *
 * TODO: update defn here, or reference to some document
 *
 * for some changepoint of interest thj and endpoints
 *
 * t_L = max(0, thj - window_size + 1)
 * t_R = min(data_count-1, thj + window_size)
 *
 * and window_size
 *
 */


double * construct_v(int data_count, int thj, int window_size, double decay_rate){
  double * v = new double[data_count];
  int tau_L = std::max(0, thj-window_size+1);
  int tau_R = std::min(data_count-1, thj+window_size);
  double const_left = -decay_rate * (pow(decay_rate, 2.0) - 1) / (pow(decay_rate, 2.0) - pow(decay_rate , (2 * (tau_L - thj) )));
  double const_right = (pow(decay_rate, 2.0) - 1) / (pow(decay_rate, (2 * (tau_R - thj))) - 1);
  for (int i = 0; i<data_count; i++){
    if(i<tau_L){
      v[i] = 0;
    }else if((i>=tau_L)&(i<=thj)){
      v[i] = const_left * pow(decay_rate, i - thj);
    }else if((i>tau_L)&(i<=tau_R)){
      v[i] = const_right * pow(decay_rate, i - thj - 1);
    }else{
      v[i] = 0;
    }
  }
  return(v);
}

// trying to test construct v
double construct_nu_norm(int data_count, int thj, int window_size, double decay_rate){
  double gam_2 = pow(decay_rate, 2.0);
  int tau_L = std::max(0, thj - window_size + 1);
  int tau_R = std::min(data_count-1, thj + window_size);
//  printf("L %i, R %i \n",tau_L,tau_R);
  double v_norm_2 = gam_2 * (gam_2 - 1)/(gam_2 - pow(decay_rate, 2 * (tau_L - thj))) + (gam_2 - 1) / (pow(decay_rate, 2 * (tau_R - thj)) - 1);
  return(v_norm_2);
}

double construct_vTy(double * y, double * v, int data_count, int thj, int window_size){
  double result = 0;
  int tau_L = std::max(0, thj-window_size+1);
  int tau_R = std::min(data_count-1, thj+window_size);
//  printf("L %i, R %i \n",tau_L,tau_R);
  for (int i = tau_L; i<=tau_R; i++){
    result += y[i]*v[i];
  }
  return(result);
}

PiecewiseSquareLoss thj_in_model(
        PiecewiseSquareLosses *cost_fwd,
        PiecewiseSquareLosses *cost_rev,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double *data_vec_rev, // flipped data
        double decay_rate, // AR1 decay parameter
        double penalty, // tuning parameter to penalize the number of spikes
        int verbose
) {

  if (window_size < 1) throw;

  int sub_f_start = std::max(thj - window_size + 1, 0); // tL
  int n_sub_f = thj - sub_f_start + 1; // number of points in forward direction
  double * sub_data_f = subset_array(data_vec, sub_f_start, thj + 1); // subsets the original data to tL:thj

  int sub_r_start = std::max(data_count - 1 - thj - window_size, 0); // tR
  int n_sub_r = data_count - 1 - thj - 1 - sub_r_start + 1; // number of points in the reverse direction
  double * sub_data_r = subset_array(data_vec_rev, sub_r_start, data_count - thj - 1); // subsets original to thj+1:tR

  double * v = construct_v(data_count, thj, window_size, decay_rate);
  double vTy = construct_vTy(data_vec, v, data_count, thj ,window_size);
  double v_norm2 = construct_vTy(v, v, data_count, thj, window_size);
  free(v); //free the vectors
  int cost_f_start = (int) std::max(thj - window_size, 0);
  int cost_r_start = (int) std::max(data_count - 1 - thj - window_size - 1, 0);

  PiecewiseSquareLoss *cost_f_piece, *cost_r_piece;
  PiecewiseSquareLoss fwd_0, rev_0;

  if (thj - window_size < 0) {
    fwd_0.piece_list.emplace_back(0, 0, 0, MACHINE_MIN, MACHINE_MAX, 0, 0); // degenerate cost section to update
    cost_f_piece = &fwd_0;
  } else {
    cost_f_piece = &(cost_fwd -> at(cost_f_start));
  }

  if (data_count - 1 - thj - window_size - 1 < 0) {
    rev_0.piece_list.emplace_back(0, 0, 0, MACHINE_MIN, MACHINE_MAX, 0, 0);
    cost_r_piece = &rev_0;
  } else {
    cost_r_piece = &(cost_rev -> at(cost_r_start));
  }


  PiecewiseBiSquareLosses fwd_2d = fpop_2d_custom_start(sub_data_f, n_sub_f,cost_f_piece, penalty, decay_rate, 1, vTy, v_norm2, verbose);
  PiecewiseBiSquareLosses rev_2d = fpop_2d_custom_start(sub_data_r, n_sub_r,cost_r_piece, penalty, decay_rate, 0, vTy, v_norm2, verbose);

  free(sub_data_f); // free pointer
  free(sub_data_r); // free pointer

  PiecewiseSquareLoss fwd_min, rev_min, c_change_at_thj;
//    printf("fwd_2d\n");
//    fwd_2d.print();
//    printf("rev_2d\n");
//    rev_2d.print();

  fwd_min = fwd_2d.min_u().get_univariate_p();
  rev_min = rev_2d.min_u().get_univariate_p();

//   printf("forward min u\n");
//    fwd_min.print();
//    printf("reverse min u\n");
//    rev_2d.min_u().print();

  c_change_at_thj.set_to_addition_of(&fwd_min, &rev_min, 0);

//  printf("change at thj \n");
//  c_change_at_thj.print();

  c_change_at_thj.add(0, 0, penalty); // add penalty at the end
  c_change_at_thj.set_prev_seg_end(1); // there is a changepoint at thj




//  printf("eval at baseline (change)= %f\n", c_change_at_thj.findCost(vTy));

  PiecewiseBiSquareLosses c_no_change_at_thj;

//  printf("forward 2d collection \n");
//  fwd_2d.print();
//
//  printf("reverse 2d collection \n");
//  rev_2d.print();

  rev_2d.rescale(0, decay_rate); //rescaling for addition - same as scaling fwd since we are minimizing over u

  //printf("reverse collection after rescaling\n");
//  printf("forward\n");
//  fwd_2d.print();
//  printf("reverse\n");
//  rev_2d.print();

  c_no_change_at_thj.set_to_addition_of(&fwd_2d, &rev_2d, 0);
  //printf("c_no_change_at_thj \n");
  //c_no_change_at_thj.print();

  PiecewiseSquareLoss c_no_change;

  c_no_change = c_no_change_at_thj.min_u().get_univariate_p();

//  printf("change at thj min over u\n");
//  c_change_at_thj.print();

  //printf("printing c_no_change...min over u\n");
  //c_no_change.print();
//  printf("eval at baseline (no change)= %f\n", c_no_change.findCost(-1.0));

//  PiecewiseBiSquareLoss c_no_change;
//  c_no_change = c_no_change_at_thj.min_u();
//  c_no_change.print();

  c_no_change.set_prev_seg_end(0); // there is no changepoint at thj

  PiecewiseSquareLoss optimal_cost_in_phi;
  optimal_cost_in_phi.set_to_min_env_of(&c_change_at_thj, &c_no_change, 0);

  //printf("final pw quadratics\n");
  //optimal_cost_in_phi.print();
  return optimal_cost_in_phi;

}

void check_selective_inference(PiecewiseSquareLoss * analytic_phi,
        int thj, // changepoint of interest
        int window_size, // size of window around thj
        int data_count, // number of data points
        double *data_vec, // original data
        double decay_rate, // decay rate of the Ca2+
        double penalty, // tuning parameter to penalize the number of spikes
        double sig, // estimated variance to limit our check range
        int verbose) {

  double v_norm2 = construct_nu_norm(data_count, thj, window_size, decay_rate);
  double * v = construct_v(data_count, thj, window_size, decay_rate);
  double vTy = construct_vTy(data_vec, v, data_count, thj, window_size);
  const double MIN = -1*std::max(10*sqrt(v_norm2*sig),fabs(vTy));
  const double MAX = std::max(10*sqrt(v_norm2*sig),fabs(vTy));


  SquareLossPieceList::iterator it;
  double phi_eval, analytic_cost, manual_cost;

  for (it = analytic_phi->piece_list.begin(); it != analytic_phi->piece_list.end(); it++) {
    phi_eval = MidMean(it -> min_mean, it -> max_mean);
    if (phi_eval > MIN && phi_eval < MAX) {
      analytic_cost = it -> getCost(phi_eval);
      // run fpop on yphi
//      printf("analytic_cost %f", analytic_cost,"\n");

      PiecewiseSquareLoss *cost_prev;
      PiecewiseSquareLosses cost_model_mat(data_count);
      PiecewiseSquareLoss start;
      start.piece_list.emplace_back(0, 0, 0, MACHINE_MIN, MACHINE_MAX, 0, 0); // degenerate cost section to update
      cost_prev = &start;

      // build cost functions for each data point, starting with start_cost cost function
      double next_data_point;
      const int fwd = 1;
      for(int data_i=0; data_i < data_count; data_i++){
          next_data_point = data_vec[data_i] - v[data_i] / v_norm2 * (vTy - phi_eval);
          cost_prev = fpop_update_i(&cost_model_mat[data_i], cost_prev, next_data_point, decay_rate, penalty, data_i, fwd, verbose);
      }
      manual_cost = cost_model_mat[data_count-1].getMinCost();
      if (ABS(manual_cost - analytic_cost) > DIFF_EPS) {
        printf("analytic cost %f, manual cost %f \n", analytic_cost, manual_cost);
        printf("analytic cost incorrect. different between analytic cost and manual cost at phi (%f)= \t %.50f", phi_eval, analytic_cost - manual_cost);
        throw std::runtime_error("analytic cost incorrect. Please report!");
      }
    }
  }
  free(v);
}

// for one-sided p-val, this computes P(X\geq vTy | X \in S) where S is the truncation set
double calc_p_value(PiecewiseSquareLoss * analytic_phi,
                    int thj, // changepoint of interest
                    int window_size, // size of window around thj
                    int data_count, // number of data points
                    double *data_vec, // original data
                    double decay_rate, // AR1 decay rate
                    double sig, // noise variance
                    bool two_sided, // whether calc 2 sided alternative or one sided (default to >)
                    double mu = 0,
                    double lower_trunc = -INFINITY) {


  double * v = construct_v(data_count, thj, window_size, decay_rate);
  double vTy = construct_vTy(data_vec, v, data_count, thj, window_size);
  double nu_norm = construct_vTy(v, v, data_count, thj, window_size);

  SquareLossPieceList::iterator it;
  free(v); // free memory
  // numerically safe
  double n1 = -INFINITY;
  double d1 = -INFINITY;
  double arg2;
  // truncation constant set to be 10||nu||sigma
  double C_stable = 0;//std::max(20*sqrt(nu_norm*sig),fabs(vTy));
  double p_val_result;
//  printf("sqrt(nu_norm * sig) %f", sqrt(nu_norm * sig));

  for (it = analytic_phi->piece_list.begin(); it != analytic_phi->piece_list.end(); it++) {
    if ((it->data_i == 1) & (it -> max_mean > lower_trunc)){ // this segment is contained and above the truncation limit

      it -> min_mean = std::max(it -> min_mean, lower_trunc);
      double a, b;
      a = pnorm_log((it -> max_mean - mu ) / sqrt(nu_norm * sig));
      b = pnorm_log((it -> min_mean - mu ) / sqrt(nu_norm * sig));
//      printf("max %f, min %f, sig %f \n", it -> max_mean, it -> min_mean,  sqrt(nu_norm * sig));


      arg2 = log_subtract(a, b);
      //printf("a %f, b %f, arg2 %f \n", a, b, arg2);

      d1 = log_sum_exp(d1, arg2);

        if (two_sided){
          if ((it->max_mean - mu) >= (ABS(vTy)-mu)) {
              arg2 = log_subtract(pnorm_log((it -> max_mean - mu) / sqrt(nu_norm * sig)),
                                  pnorm_log(std::max(it -> min_mean - mu, ABS(vTy)-mu) / sqrt(nu_norm * sig)));
              n1 = log_sum_exp(n1, arg2);
          }
          if ((it->min_mean - mu) <= (-1 * ABS(vTy)-mu)) {
              arg2 = log_subtract(pnorm_log(std::min(it -> max_mean - mu , -ABS(vTy) - mu) / sqrt(nu_norm * sig)),
                                  pnorm_log((it -> min_mean - mu) / sqrt(nu_norm * sig)));
              n1 = log_sum_exp(n1, arg2);
          }
      } else { //one sided p-value
          if ((it->max_mean-mu) >= (vTy-mu)) {
              arg2 = log_subtract(pnorm_log((it -> max_mean - mu) / sqrt(nu_norm * sig)),
                                  pnorm_log(std::max(it -> min_mean - mu, (vTy-mu)) / sqrt(nu_norm * sig)));
              n1 = log_sum_exp(n1, arg2);
          }
      }

    }
  }


//   printf("n1 %f d1 %f\n",n1,d1);

    double result = exp(n1-d1);

    if (result!=result){
          p_val_result = 0.0; // numerical instability; long tail
      }else{
          p_val_result = exp(n1 - d1);
      }
//    printf("p_val_result %f\n",p_val_result);
    return (p_val_result);
}


double calc_surv_prob(PiecewiseSquareLoss * analytic_phi,
                    double mu, // mean of the truncated gaussian
                    double vTy, //
                    double nu_norm,
                    double sig,
                    double lower_trunc) {

    SquareLossPieceList::iterator it;
    double p_val_result;
    // numerically safe
    double n1 = -INFINITY;
    double d1 = -INFINITY;
    double arg2;
    for (it = analytic_phi->piece_list.begin(); it != analytic_phi->piece_list.end(); it++) {
        if ((it->data_i == 1)  & (it -> max_mean > lower_trunc)) { // this segment is contained
            it -> min_mean = std::max(it -> min_mean, lower_trunc);
            double a, b;
            a = pnorm_log((it -> max_mean - mu) / sqrt(nu_norm * sig));
            b = pnorm_log((it -> min_mean - mu) / sqrt(nu_norm * sig));
            arg2 = log_subtract(a, b);
            d1 = log_sum_exp(d1, arg2);
            if ((it->max_mean - mu) >= (vTy-mu)) {
                arg2 = log_subtract(pnorm_log((it -> max_mean - mu ) / sqrt(nu_norm * sig)),
                                    pnorm_log(std::max(it -> min_mean - mu, vTy - mu) / sqrt(nu_norm * sig)));
                n1 = log_sum_exp(n1, arg2);
            }
        }
    }

    double result = exp(n1-d1);

    if (result!=result){
        p_val_result=0.0;
    }else{
        p_val_result = exp(n1 - d1);
    }
    return (exp(n1 - d1));

}


double tn_lower_surv(PiecewiseSquareLoss * analytic_phi,
                     double mu, // mean of the truncated gaussian
                     double vTy,
                     double nu_norm,
                     double sig, // noise variance
                     double alpha_1,
                     double lower_trunc){

     double result = calc_surv_prob(analytic_phi,
             mu,
             vTy,
             nu_norm,
             sig,
             lower_trunc)-(alpha_1);

     return(result);
}


double tn_upper_surv(PiecewiseSquareLoss * analytic_phi,
                     double mu, // mean of the truncated gaussian
                     double vTy,
                     double nu_norm,
                     double sig, // noise variance
                     double alpha_2,
                     double lower_trunc){
    double result = calc_surv_prob(analytic_phi,
            mu,
            vTy,
            nu_norm,
            sig,
            lower_trunc)-(1-alpha_2);
    return(result);

}



// compute the confidence interval for the mean parameter of a normal truncated to a union of intervals
std::vector<double> compute_CI(PiecewiseSquareLoss * analytic_phi,
                    int thj, // changepoint of interest
                    int window_size, // size of window around thj
                    int data_count, // number of data points
                    double *data_vec, // original data
                    double decay_rate, // AR1 decay rate
                    double sig, // noise variance
                    double alpha = 0.05, // desired type I error rate control, default to 0.05
                    double lower_trunc = -INFINITY// lower_trunc to guide numerically stable CI computation
                   ){

    //BISECTION_EPS
    double alpha_1 = alpha/2.0;
    double alpha_2 = alpha-alpha_1;
    double lower_CI, upper_CI;
    double xL_1 = 0;
    double xR_1 = 1;
    double xL_2 = 0;
    double xR_2 = 1;
    double f_lower_neg, f_upper_neg, delta, f_lower_pos, f_upper_pos;

    // recycle these computation
    double * v = construct_v(data_count, thj, window_size, decay_rate);
    double vTy = construct_vTy(data_vec, v, data_count, thj, window_size);
    double nu_norm = construct_vTy(v, v, data_count, thj, window_size);
    free(v); // free memory


    f_lower_neg = tn_lower_surv(analytic_phi,
            xL_1,
            vTy,
            nu_norm,
            sig, // noise variance
            alpha_1,
            lower_trunc); // ideally this is <= 0

    if (f_lower_neg < 0){
        f_lower_pos = tn_lower_surv(analytic_phi,
                                    xL_2,
                                    vTy,
                                    nu_norm,
                                    sig, // noise variance
                                    alpha_1,
                                    lower_trunc);
        while (f_lower_pos < 0 && xL_2 <= 1E3){
            delta = 0.01*std::max(1e-3, fabs(f_lower_pos));
            xL_2 = xL_2+delta;
            f_lower_pos = tn_lower_surv(analytic_phi,
                                        xL_2,
                                        vTy,
                                        nu_norm,
                                        sig, // noise variance
                                        alpha_1,
                                        lower_trunc);
        }
    }else{
        while (f_lower_neg > 0 && xL_1 >= -1E3){

            delta = 0.01*std::max(1e-3, fabs(f_lower_neg));
            xL_1 = xL_1-delta;
            f_lower_neg = tn_lower_surv(analytic_phi,
                                    xL_1,
                                    vTy,
                                    nu_norm,
                                    sig, // noise variance
                                    alpha_1,
                                    lower_trunc);
        }
    }

    f_upper_neg = tn_upper_surv(analytic_phi,
                            xR_1,
                            vTy,
                            nu_norm,
                            sig, // noise variance
                            alpha_2,
                            lower_trunc);

    if (f_upper_neg < 0){
        f_upper_pos = tn_upper_surv(analytic_phi,
                                    xR_2,
                                    vTy,
                                    nu_norm,
                                    sig, // noise variance
                                    alpha_2,
                                    lower_trunc);

        while (f_upper_pos < 0 && xR_2 <= 1E3){
            delta = 0.01*std::max(1e-3, (double) fabs(f_upper_pos));
            xR_2 = xR_2+delta;
            f_upper_pos = tn_upper_surv(analytic_phi,
                                        xR_2,
                                        vTy,
                                        nu_norm,
                                        sig, // noise variance
                                        alpha_2,
                                        lower_trunc);
        }

    }else{
        while (f_upper_neg>0 && xR_1 >= -1E3){
            delta = 0.01*std::max(1e-3, (double) fabs(f_upper_neg));
            xR_1 = xR_1-delta;
            f_upper_neg = tn_upper_surv(analytic_phi,
                                        xR_1,
                                        vTy,
                                        nu_norm,
                                        sig, // noise variance
                                        alpha_2,
                                        lower_trunc);
        }
    }

    // xL_1, xL_2, xR_1, xR_2 are ready to pass into the root finding eqs
    // Bisection to find lower_CI
//    printf("bisection %f, %f, %f, %f \n",xL_1,xL_2,xR_1,xR_2);
    double f_midpoint, f_1;
    while (xL_2 - xL_1 >= BISECTION_EPS){
        double mid_point = (xL_1 + xL_2)/2.0;
        f_1 = tn_lower_surv(analytic_phi,
                            xL_1,
                            vTy,
                            nu_norm,
                            sig, // noise variance
                            alpha_1,
                            lower_trunc);

        f_midpoint = tn_lower_surv(analytic_phi,
                                          mid_point,
                                          vTy,
                                          nu_norm,
                                          sig,
                                          alpha_1,
                                   lower_trunc);

        if (f_1*f_midpoint>=0){ // same sign
            xL_1 = mid_point;
        }else{
            xL_2 = mid_point;
        }
    }

    lower_CI = (xL_1 + xL_2)/2.0;

    // Bisection to find upper_CI
    while (xR_2 - xR_1 >= BISECTION_EPS){
        double mid_point = (xR_1 + xR_2)/2.0;
        f_1 = tn_upper_surv(analytic_phi,
                            xR_1,
                            vTy,
                            nu_norm,
                            sig,
                            alpha_2,
                            lower_trunc);

        f_midpoint = tn_upper_surv(analytic_phi,
                                   mid_point,
                                   vTy,
                                   nu_norm,
                                   sig, // noise variance
                                   alpha_2,
                                   lower_trunc);

        if (f_1*f_midpoint>=0){ // same sign
            xR_1 = mid_point;
        }else{
            xR_2 = mid_point;
        }
    }

    upper_CI = (xR_1 + xR_2)/2.0;

    std::vector<double> CI_result = {lower_CI,upper_CI};
//    printf("lower_CI %f upper_CI %f \n",lower_CI,upper_CI);
    return (CI_result);
}