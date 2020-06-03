#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "fpop.h"
#include "utils.h"
#include "selective_inference.h"
#include "fpop_inference.h"
#include "test_utils.h" // for unit test purposes

typedef std::vector<double> VecDouble;
typedef std::vector<int> VecInt;

VecDouble read_data_vec_double(const std::string filename, const int data_count) {
  const int max_line_size = 100;
  char in[max_line_size];

  std::ifstream fsnumbers;
  fsnumbers.open(filename);

  VecDouble data_vec;
  for (int i = 0; i < data_count; i++) {
    fsnumbers >> in;
    data_vec.emplace_back(std::stod(in));
  }
  return data_vec;
}



const static double EPS = 1e-6; // numerical tolerance
int tests_total = 0;
int tests_fail = 0;

bool vector_within_eps(double * vec_a, double * vec_b, double eps) {

    int a_length = sizeof( vec_a )/sizeof( vec_a[0] );
    int b_length = sizeof( vec_b )/sizeof( vec_b[0] );

    assert(a_length == b_length);

    for (int i = 0; i < a_length; i++) {
        double curr_diff = std::abs(vec_a[i] - vec_b[i]);
        if (curr_diff > eps){
            return false;
        }
    }
    return true;
}

bool double_within_eps (double a, double b, double eps) {
    if (std::abs(a - b) > eps) {
        return false;
    }
    return true;
}


double naive_nuTy (double * vec_a, double * vec_b, int data_count) {
    double result = 0;
    for (int i = 0; i < data_count; i++) {
        result += vec_a[i]*vec_b[i];
    }
    return(result);
}


void toy_example_1(){

    const int data_count = 6;
    double y[data_count] =  {8, 4, 2, 5, 2.5, 1.25};

    // print out y data
    double penalty = 1;
    for (auto x : y) {printf("%f \t", x);}
    printf("\n");

    int thj = 2;
    int window_size = 1;
    double decay_rate = 0.5;
    const double sig = 1;

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, -INFINITY, INFINITY);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, -INFINITY, INFINITY);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate, penalty, thj, window_size, sig);

    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);

    printf("Cost model, short data, 1 true cp \n");
    out.model.print();
    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

void toy_example_2(){

    const int data_count = 11; // # of data points
//  const std::string filename = "/Users/jewellsean/Desktop/test.csv";
//  VecDouble y = read_data_vec_double(filename, data_count);
    double y[data_count] = {8, 4, 2, 1, 0.5, 8, 4, 2, 1, 0.5, 10};

    // print out y data
    double penalty = 0.5;
    for (auto x : y) {printf("%f \t", x);}
    printf("\n");

    int thj = 4;
    int window_size = 3;
    double decay_rate = 0.5;
    const double sig = 1;

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, -INFINITY, INFINITY);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, -INFINITY, INFINITY);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate, penalty, thj, window_size, sig);

    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);

    printf("Cost model, longer data, 2 true cp \n");
    out.model.print();
    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

void toy_example_3(){

    const int data_count = 11; // # of data points
//  const std::string filename = "/Users/jewellsean/Desktop/test.csv";
//  VecDouble y = read_data_vec_double(filename, data_count);
    double y[data_count] = {8, 4, 2, 1, 0.5, 8, 4, 2, 1, 0.5, 10};

    // print out y data
    double penalty = 0.5;
    for (auto x : y) {printf("%f \t", x);}
    printf("\n");

    int thj = 10;
    int window_size = 1;
    double decay_rate = 0.5;
    const double sig = 1;

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, -INFINITY, INFINITY);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, -INFINITY, INFINITY);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate, penalty, thj, window_size, sig);

    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);

    printf("Cost model, longer data, 2 true cp \n");
    out.model.print();
    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

int main(int argc, char *argv[]) {
//  const std::string filename = "/Users/jewellsean/Desktop/test.csv";
//  VecDouble y = read_data_vec_double(filename, data_count);

 toy_example_1();
// toy_example_2();
// toy_example_3();
 return 0;

}



// CODE for printing out v
//  double * v_test = construct_v(data_count, thj, window_size, decay_rate);

//
//  for(int data_i=0; data_i < data_count; data_i++){
//    printf("%f \t", v_test[data_i]);
//    printf("\n");
//  }
//
//  double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
//  double v_norm2 = construct_nu_norm(data_count, thj, window_size, decay_rate);
//  printf("vTy %f\n",vTy);
//  printf("vTv %f\n",v_norm2);




int main2(int argc, char *argv[]) {
  const int data_count = 11; // # of data points
//  const std::string filename = "/Users/jewellsean/Desktop/test.csv";
//  VecDouble y = read_data_vec_double(filename, data_count);

  double y[data_count] = {8, 4, 2, 1, 0.5, 8, 4, 2, 1, 0.5, 10};

//  printf("Array length %i \n",sizeof( y )/sizeof( y[0] ));

  double penalty = 1;
//  for (auto x : y) {printf("%f \t", x);}
//  printf("\n");

  int thj = 3;
  int window_size = 2;
  double decay_rate = 0.5;

  // the correct v constructed using R
  double v_reference[data_count] = {0.0, 0.0,  0.0, -0.2, -0.1,  0.8,  0.4,  0.0,  0.0,  0.0,  0.0};
  double * v_test = construct_v(data_count, thj, window_size, decay_rate);

  // TEST 0: correctness of contrast vector
  TEST(vector_within_eps(v_reference, v_test, EPS));
  // TEST 1: correctness of contrast vector: should return close 0 with no cp
  double * null_v_test = construct_v(data_count, 2, window_size, decay_rate);
  double null_vTy = construct_vTy(y, null_v_test, data_count, 2 ,window_size);
  TEST(double_within_eps(null_vTy, 0.0, EPS));

  // TEST 2: test the summation vTy
  double test_vTy = naive_nuTy(y, v_test, data_count);
  double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
  TEST(double_within_eps(test_vTy, vTy, EPS));

  // TEST 3: v_norm2
  double test_vTv = naive_nuTy(v_test, v_test, data_count);
  double v_norm_2 = construct_nu_norm(data_count, thj, window_size, decay_rate);
  TEST(double_within_eps(test_vTv, v_norm_2, EPS));

//  for (int i = 0; i < data_count; i++){
//    printf("%f \t", null_v_test[i]);
//  }
//
//  printf("\n vTy = %f\n", null_vTy);

  return 0;
}
