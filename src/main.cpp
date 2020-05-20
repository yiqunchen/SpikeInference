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


int main(int argc, char *argv[]) {
  const int data_count = 11; // # of data points
//  const std::string filename = "/Users/jewellsean/Desktop/test.csv";
//  VecDouble y = read_data_vec_double(filename, data_count);

  double y[data_count] = {8, 4, 2, 1, 0.5, 8, 4, 2, 1, 0.5, 10};

  double penalty = 1;
  for (auto x : y) {printf("%f \t", x);}
  printf("\n");

  int thj = 5;
  int window_size = 2;
  double decay_rate = 0.5;

  double * v_test = construct_v(data_count, thj, window_size, decay_rate);
  double vTy = construct_vTy(y, v_test, data_count, thj ,window_size);

  for (int i = 0; i < data_count; i++){
    printf("%f \t", v_test[i]);
  }

  printf("\n vTy = %f\n", vTy);

//  PiecewiseSquareLosses out = fpop(y, data_count, gam, penalty, 0.0, INFINITY);
//
//  double * data_vec_rev = reverse_data(y, data_count);
//  PiecewiseSquareLosses cost_model_rev = fpop(data_vec_rev, data_count, 1.0 / gam, penalty, 0.0, INFINITY);
//
//
////  for( auto x : out) {x.print();}
//  printf("Estimated changepoints in forward direction\n");
//  std::list<int> change_pts = extract_changepoints(out, data_count);
//  for (auto x : change_pts) {printf("%d\t", x);}
//
//  printf("\n Estimated changepoints in reverse direction\n");
//  std::list<int> change_pts_rev = extract_changepoints(cost_model_rev, data_count);
//  for (auto x : change_pts_rev) {printf("%d\t", x);}

}

