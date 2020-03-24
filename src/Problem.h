/* 
  Problem.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#include <stdlib.h>
#include<stdio.h>
#include <random>
#include<iostream>
#include <sstream>

#include <iomanip>
#include <cassert>
#include <string.h>
#include <vector>
#include<math.h>
#include <sys/stat.h>
#include <map>

using std::cout;
using std::endl;
using std::string;

class Problem {
 public:
  Problem(string test_func, int dim, string variable_order, string func_type, int bbob_func_id);
  ~Problem();
  
  double eval(double *x);
  void setInstanceData();
  double errorValBSF();
  bool updateBSFsolution(double f, double *x);
  double sphere(double *x);
  double rosenbrock(double *x);
  double rastrigin(double *x);
  string int2string(int i);
  int string2int(const string& str);
  double string2double(const string& str);
    
  inline int getDimSize(){return this->dim;}
  inline double getUpperBound(){return this->upper_bound;}
  inline double getLowerBound(){return this->lower_bound;}
  inline string getTestFuncName(){return this->test_func;}
  
 private:
  string test_func;
  int dim;
  string variable_order;       
  string func_type;
  double val2rearch;
  double bsf_obj_val;
  double *bsf_solution;  
  double opt_obj_val;
  double *shifted_opt_pos;
  int *random_permutation;
  double ***d2_rot_matrix;
  double *tmp_vec;
  double *translated_x;
  double **rot_matrix;
  double upper_bound;
  double lower_bound;  
};

#endif
