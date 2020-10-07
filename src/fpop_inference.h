#include "funPieceList.h"
#include <stdexcept>

class FpopInference {
public:
    double pval;
    double approximation_error;
    PiecewiseSquareLoss model;
    int thj;
    std::vector<double> confidence_interval;
    FpopInference(double p, double a, PiecewiseSquareLoss m, int t, std::vector<double> confidence_interval);
    FpopInference(double p, double a, PiecewiseSquareLoss m, int t);
};

FpopInference fpop_analytic_inference_recycle(PiecewiseSquareLosses * cost_model_fwd,
                                       PiecewiseSquareLosses * cost_model_rev,
                                       double * data_vec, int data_count,
                                       double decay_rate,
                                       double penalty,
                                       int thj,
                                       int window_size,
                                       double sig,
                                       bool return_ci,
//                                       bool two_sided,
//                                       double alpha,
//                                       double mu);
                                       double alpha);

