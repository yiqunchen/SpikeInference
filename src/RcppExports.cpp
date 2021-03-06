// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fpop_interface2
List fpop_interface2(NumericVector data, double gam, double penalty, double min_mean, double max_mean, NumericVector cost_mat_r, IntegerVector end_vec_r, NumericVector mean_vec_r, IntegerVector intervals_mat_r);
RcppExport SEXP _SpikeInference_fpop_interface2(SEXP dataSEXP, SEXP gamSEXP, SEXP penaltySEXP, SEXP min_meanSEXP, SEXP max_meanSEXP, SEXP cost_mat_rSEXP, SEXP end_vec_rSEXP, SEXP mean_vec_rSEXP, SEXP intervals_mat_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< double >::type min_mean(min_meanSEXP);
    Rcpp::traits::input_parameter< double >::type max_mean(max_meanSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cost_mat_r(cost_mat_rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type end_vec_r(end_vec_rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean_vec_r(mean_vec_rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type intervals_mat_r(intervals_mat_rSEXP);
    rcpp_result_gen = Rcpp::wrap(fpop_interface2(data, gam, penalty, min_mean, max_mean, cost_mat_r, end_vec_r, mean_vec_r, intervals_mat_r));
    return rcpp_result_gen;
END_RCPP
}
// fpop_inference_interface_recycle
List fpop_inference_interface_recycle(NumericVector data, double decay_rate, double penalty, int window_size, double sig, int return_dev, bool return_ci, bool two_sided, double alpha, double mu, double lower_trunc);
RcppExport SEXP _SpikeInference_fpop_inference_interface_recycle(SEXP dataSEXP, SEXP decay_rateSEXP, SEXP penaltySEXP, SEXP window_sizeSEXP, SEXP sigSEXP, SEXP return_devSEXP, SEXP return_ciSEXP, SEXP two_sidedSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP lower_truncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type decay_rate(decay_rateSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< int >::type return_dev(return_devSEXP);
    Rcpp::traits::input_parameter< bool >::type return_ci(return_ciSEXP);
    Rcpp::traits::input_parameter< bool >::type two_sided(two_sidedSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lower_trunc(lower_truncSEXP);
    rcpp_result_gen = Rcpp::wrap(fpop_inference_interface_recycle(data, decay_rate, penalty, window_size, sig, return_dev, return_ci, two_sided, alpha, mu, lower_trunc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpikeInference_fpop_interface2", (DL_FUNC) &_SpikeInference_fpop_interface2, 9},
    {"_SpikeInference_fpop_inference_interface_recycle", (DL_FUNC) &_SpikeInference_fpop_inference_interface_recycle, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpikeInference(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
