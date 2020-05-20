#include "funPieceList.h"
#include <stdexcept>

PiecewiseSquareLosses fpop(double *, int, double, double, double, double);

PiecewiseSquareLosses fpop_custom(double *, int, PiecewiseSquareLoss *, double, double, bool, int);

void decode_fpop(PiecewiseSquareLosses, int, double *, int *, double *, int *, double gam);

void check_min_operation(PiecewiseSquareLoss *, PiecewiseSquareLoss *, PiecewiseSquareLoss *, double, int);

PiecewiseSquareLoss * fpop_update_i(PiecewiseSquareLoss *, PiecewiseSquareLoss *, double, double, double, int, bool, int);

PiecewiseBiSquareLosses fpop_2d_custom_start(double *data_vec, int data_count,
                                             PiecewiseSquareLoss *start_cost,
                                             double penalty,
                                             double decay_rate,
                                             bool forward,
                                             double nuTy,
                                             double norm_constant,
                                             int verbose);

std::list<int> extract_changepoints(PiecewiseSquareLosses cost_model_mat, int data_count);