#include <vector>
#include <stdio.h>
#include <stdlib.h>    /* srand, rand */
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
#include <cassert>

typedef std::vector<double> VecDouble;
typedef std::vector<int> VecInt;

//VecDouble read_data_vec_double(const std::string filename, const int data_count) {
//  const int max_line_size = 20;
//  std::string in;
//
//  std::ifstream fsnumbers;
//  fsnumbers.open(filename);
//
////    fsnumbers.getline(in, max_line_size);
//  VecDouble data_vec;
//  //for (int i = 0; i < data_count; i++) {
//  while (std::getline(fsnumbers, in)){
//      data_vec.emplace_back(std::stod(in));
//  }
//  return data_vec;
//}

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
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate, penalty, thj,
            window_size, sig, true, false, 0.05, 0, -INFINITY);
    printf("true vTc %f \n", y[thj+1]-y[thj]);


//    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
//    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
//
//    printf("Cost model, longer data, 2 true cp \n");
//    out.model.print();
//    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

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

    int thj = 9;
    int window_size = 2;
    double decay_rate = 0.5;
    const double sig = 1;

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate,
            penalty, thj, window_size, sig, true, false, 0.05, 0, -INFINITY);

//    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
//    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
//
//    printf("Cost model, longer data, 2 true cp \n");
//    out.model.print();
//    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

    printf("true vTc %f \n", y[thj+1]-y[thj]);


}

void toy_example_3(){

  const int data_count = 4; // # of data points
  double y[data_count] = {8, 4, 2, 5};

  // print out y data
  double penalty = 1;
  for (auto x : y) {printf("%f \t", x);}
  printf("\n");

  int thj = 2;
  int window_size = 2;
  double decay_rate = 0.5;
  const double sig = 1;

  // forward pass
  PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
  // backward pass
  double *data_vec_rev = reverse_data(y, data_count);
  PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
  FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate,
          penalty, thj, window_size, sig, true, false, 0.05, 0, -INFINITY);

  printf("true vTc %f \n", y[thj+1]-y[thj]);

    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
//  for(int data_i=0; data_i < data_count; data_i++){
//    printf("%f \n", v_test[data_i]);
//  }


  double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
  double v_norm2 = construct_vTy(v_test, v_test, data_count, thj, window_size);
  double v_norm2_incorrect = construct_nu_norm(data_count, thj, window_size, decay_rate);
//  printf("vTy %f\n",vTy);
//  printf("vTv %f\n",v_norm2);
//  printf("vTv false %f\n",v_norm2_incorrect);


//  double phi_eval = 0;
//  for(int data_i=0; data_i < data_count; data_i++){
//      double next_data_point = y[data_i] - v_test[data_i] / v_norm2 * (vTy - phi_eval);
//      printf("%f \t", next_data_point);
//    }

//
//
//  printf("Cost model, longer data, 1 true cp \n");
//  out.model.print();
//  printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

void toy_example_4(){

    const int data_count = 5; // # of data points
    double y[data_count] = {8, 4, 6, 3, 1.5};

    // print out y data
    double penalty = 1;
    for (auto x : y) {printf("%f \t", x);}
    printf("\n");

    int thj = 2;
    int window_size = 2;
    double decay_rate = 0.5;
    const double sig = 1;


    double * v_test = construct_v(data_count, thj, window_size, decay_rate);

    for(int data_i=0; data_i < data_count; data_i++){
        printf("%f \t", v_test[data_i]);
    }
    printf("\n");

    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
    double vTv = construct_vTy(v_test, v_test, data_count, thj, window_size);
    printf("vTy %f vTv %f\n", vTy,vTv);

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate,
            penalty, thj, window_size, sig, true, false, 0.05, 0, -INFINITY);


    printf("true vTc %f \n", y[thj+1]-y[thj]);




//
//    printf("Cost model, longer data, 1 true cp \n");
//    out.model.print();
//    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

void toy_example_5(){

    const int data_count = 6; // # of data points
    double y[data_count] = {3, 2.7, 2.43, 2.18, 2.7, 2.43};

    // print out y data
    double penalty = 0.1;
    for (auto x : y) {printf("%f \t", x);}
    printf("\n");

    int thj = 3;
    int window_size = 2;
    double decay_rate = 0.9;
    const double sig = 0.5;

    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate,
            penalty, thj, window_size, sig, true, false, 0.05, 0, -INFINITY);
    printf("p val %f", out.pval);
//    double * v_test = construct_v(data_count, thj, window_size, decay_rate);

//    for(int data_i=0; data_i < data_count; data_i++){
//        printf("%f \t", v_test[data_i]);
//        printf("\n");
//    }
//    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);

//    printf("true vTc %f \n", y[thj+1]-decay_rate*y[thj]);


//    printf("Cost model, longer data, 1 true cp \n");
//    out.model.print();
//    printf("Cost on original data = \t %f\n", out.model.findCost(vTy));

}

// random test
void random_example_test(int T, double decay_rate, double spike_rate, float sigma,
        bool noisy, double penalty, int window_size, int random_seed,
        bool parallel) {

    assert(decay_rate < 1 & decay_rate > 0);
    assert(spike_rate > 0);
    assert(T > 0);
    assert(sigma > 0);


    std::default_random_engine generator(random_seed);
    std::poisson_distribution<int> poisson(spike_rate);
    std::normal_distribution<double> normal(0.0, sigma);

    double y[T];
    double c[T];
    y[0] = 1.0;
    c[0] = 1.0;
    int thj;


    //double vTv = construct_nu_norm(T, thj, window_size, decay_rate);

    //printf("vTv %f \n", vTv);


    for (int data_i = 1; data_i < T; data_i++) {
        int spike_T = poisson(generator);
        c[data_i] = decay_rate * c[data_i - 1] + spike_T;

        if (noisy) {
            double noise_T = normal(generator);
            y[data_i] = c[data_i]+noise_T;
        }else{
            y[data_i] = c[data_i];
        }
    }

    PiecewiseSquareLosses cost_model_fwd = fpop(y, T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    std::list<int> ll = extract_changepoints(cost_model_fwd, T);
    ll.pop_front(); // we don't want to test the first loc at -1
    // convert the cp into vector so we could use openmp
    std::vector<int> ll_vec(ll.begin(), ll.end());
    std::list<int>::iterator j;
    int count = 0;
    double *data_vec_rev = reverse_data(y, T);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    free(data_vec_rev); //free memory
//    if (parallel){
//        omp_set_num_threads(4); //make this an argument?
//        #pragma omp parallel for
//        for (int k = 0; k < ll_vec.size(); ++k) {
//            count += 1;
//            thj = ll_vec[k]; // get random spike location to test
//            printf("currently testing %i th location at %i\n", count, thj);
//            // forward pass
//            PiecewiseSquareLosses cost_model_fwd = fpop(y, T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
//            // backward pass
//            double *data_vec_rev = reverse_data(y, T);
//            PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
//            FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, T, decay_rate, penalty,
//                                                                thj, window_size, sigma * sigma);
//
//        }
//    }else{
    double total_test = 0, cover_test = 0;
    for (j = ll.begin(); j != ll.end(); ++j) {
        count += 1;
        thj = *j; // get estimated spike location to test
        //printf("currently testing %i th location at %i\n", count, thj);
        // forward pass
        // backward pass
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, T, decay_rate, penalty,
                                                            thj, window_size, sigma * sigma, true, false,
                                                            0.05, 0, 0);

        double * v = construct_v(T, thj, window_size, decay_rate);
        double vTc = construct_vTy(v, c, T, thj, window_size);
        free(v); //
        //printf("true spike jump %f \n", vTc);
        total_test += 1;
        if (out.confidence_interval[0]<=vTc & out.confidence_interval[1]>=vTc){
            cover_test += 1;
        }
    }
    printf("total %f cover %f rate %f\n", total_test, cover_test, cover_test/total_test);
//}
}


void specific_example_1() {

    float decay_rate = 0.95;
    float penalty = 0.3;
    int T = 100000;
    int thj;
    int window_size = 2;
    float sigma = 0.15;
    vector<double> y_example = read_data_vec_double("/Users/tonyyiqunchen/Desktop/test.csv",
                                                    T);

    PiecewiseSquareLosses cost_model_fwd = fpop(&y_example[0], T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    std::list<int> ll = extract_changepoints(cost_model_fwd, T);
    ll.pop_front(); // we don't want to test the first loc at -1
    // convert the cp into vector so we could use openmp
    std::vector<int> ll_vec(ll.begin(), ll.end());
    std::list<int>::iterator j;
    int count = 0;

    // backward pass
    double *data_vec_rev = reverse_data(&y_example[0], T);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    for (j = ll.begin(); j != ll.end(); ++j) {
        count += 1;
        thj = *j; // get estimated spike location to test
        printf("currently testing %i th location at %i\n", count, thj);
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, &y_example[0], T, decay_rate, penalty,
                                                            thj, window_size, sigma * sigma, true, false, 0.05, 0, -INFINITY);

    }
}


void specific_example_2() {

    float decay_rate = 0.95;
    float penalty = 0.03;
    int T = 10000;
    int thj;
    int window_size = 4;
    float sigma = 0.15;
    vector<double> y_example = read_data_vec_double("/Users/tonyyiqunchen/Desktop/calcium_test_2.csv",
                                                    T);

    PiecewiseSquareLosses cost_model_fwd = fpop(&y_example[0], T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    std::list<int> ll = extract_changepoints(cost_model_fwd, T);
    ll.pop_front(); // we don't want to test the first loc at -1
    // convert the cp into vector so we could use openmp
    std::vector<int> ll_vec(ll.begin(), ll.end());
    std::list<int>::iterator j;
    int count = 0;

    // backward pass
    double *data_vec_rev = reverse_data(&y_example[0], T);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    for (j = ll.begin(); j != ll.end(); ++j) {
        count += 1;
        thj = *j; // get estimated spike location to test
        printf("currently testing %i th location at %i\n", count, thj);
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, &y_example[0], T,
                                                                decay_rate, penalty,thj, window_size, sigma * sigma,
                                                                true, false, 0.05, 0, -INFINITY);

    }
}


void specific_example_3() {

    float decay_rate = 0.9857143;//0.9857143;
    float penalty = 15;//63.09573;
    int T = 31962;
    int thj;
    int window_size = 5;
    float sigma = sqrt(0.0006845397);

    vector<double> y_example = read_data_vec_double("/Users/tonyyiqunchen/Desktop/Ca_nan_example.csv",T);
   // vector<double> y_example = read_data_vec_double("/Users/tonyyiqunchen/Desktop/calcium_ci_test.csv",T);

    PiecewiseSquareLosses cost_model_fwd = fpop(&y_example[0], T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    std::list<int> ll = extract_changepoints(cost_model_fwd, T);
    ll.pop_front(); // we don't want to test the first loc at -1
    // convert the cp into vector so we could use openmp
    std::vector<int> ll_vec(ll.begin(), ll.end());
    std::list<int>::iterator j;
    int count = 0;

    // backward pass
    double *data_vec_rev = reverse_data(&y_example[0], T);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    free(data_vec_rev);

    printf("ll begin %f", ll.begin());
    for (j = ll.begin(); j != ll.end(); ++j) {
        count += 1;
        thj = *j; // get estimated spike location to test
        printf("currently testing %i th location at %i\n", count, thj);
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, &y_example[0], T,
                                                            decay_rate, penalty, thj, window_size, sigma * sigma,
                                                            true, true, 0.05, 0, -INFINITY);
        printf("p value %f \n", out.pval);

    }
}




void specific_example_4() {

    float decay_rate = 0.98;
    float penalty = 3.826904;
    int T = 10000;
    int thj;
    int window_size = 10;
    float sigma = 1;
    vector<double> y_example = read_data_vec_double("/Users/tonyyiqunchen/Desktop/calcium_debug_10_06.csv",
                                                    T);

    PiecewiseSquareLosses cost_model_fwd = fpop(&y_example[0], T, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    std::list<int> ll = extract_changepoints(cost_model_fwd, T);
    ll.pop_front(); // we don't want to test the first loc at -1
    // convert the cp into vector so we could use openmp
    std::vector<int> ll_vec(ll.begin(), ll.end());
    std::list<int>::iterator j;
    int count = 0;

    // backward pass
    double *data_vec_rev = reverse_data(&y_example[0], T);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, T, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    for (j = ll.begin(); j != ll.end(); ++j) {
        count += 1;
        thj = *j; // get estimated spike location to test
        printf("currently testing %i th location at %i\n", count, thj);
        FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, &y_example[0], T,
                                                            decay_rate, penalty, thj, window_size,
                                                            sigma * sigma, true, false, 0.05, 0, 0);
        printf("p value %f \n", out.pval);

    }
}



void paper_example(){

    const int data_count = 4; // # of data points
    double y[data_count] = {8, 4, 2, 3};

    // print out y data
    double penalty = 0.1;
    //for (auto x : y) {printf("%f \t", x);}
    //printf("\n");

    int thj = 1;
    int window_size = 1;
    double decay_rate = 0.5;
    const double sig = 1;

    double * v_test = construct_v(data_count, thj, window_size, decay_rate);

    for(int data_i=0; data_i < data_count; data_i++){
        printf("%f \t", v_test[data_i]);
    }
    printf("\n");

//    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
//    double vTv = construct_vTy(v_test, v_test, data_count, thj, window_size);
//    printf("vTy %f vTv %f\n", vTy,vTv);


    // forward pass
    PiecewiseSquareLosses cost_model_fwd = fpop(y, data_count, decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);
    // backward pass
    double *data_vec_rev = reverse_data(y, data_count);
    PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1 / decay_rate, penalty, MACHINE_MIN, MACHINE_MAX);

    FpopInference out = fpop_analytic_inference_recycle(&cost_model_fwd, &cost_model_rev, y, data_count, decay_rate, penalty,
            thj, window_size, sig, true, false, 0.05, 0, 0);
    printf("p value %f \n", out.pval);
    printf("CI %f, %f \n", out.confidence_interval[0], out.confidence_interval[1]);
//    printf("true vTc %f \n", y[thj+1]-y[thj]);

//    double * v_test = construct_v(data_count, thj, window_size, decay_rate);
//    for(int data_i=0; data_i < data_count; data_i++){
//    printf("%f \n", v_test[data_i]);
//  }

//    double vTy = construct_vTy(y, v_test, data_count, thj, window_size);
//    double v_norm2 = construct_vTy(v_test, v_test, data_count, thj, window_size);
    //printf("%f %f \n",vTy, v_norm2);
}

int main(int argc, char *argv[]) {
 //toy_example_1();
 //toy_example_2();
 //toy_example_3();
 //toy_example_4();
 toy_example_5();
//printf("exp %f \n", exp(-1*50));
 int test_times = 1;
 //paper_example();
 //for (int i = 0; i < test_times; i++){
 //    random_example_test(100, 0.95, 0.01, 1,true, 1, 20, i, false);
 // }
 //specific_example_4();
//specific_example_4();
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

//
//  TODO: We need a better way to test the correctness of the contrast vector
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

  return 0;
}
