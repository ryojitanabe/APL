/* 
  DE.h
  Author: Ryoji Tanabe 2020
 */

#ifndef _DE_H_
#define _DE_H_

#include "Problem.h"

#include <stdlib.h>
#include<stdio.h>
#include <random>
#include<iostream>
#include <sstream>
#include<fstream>
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
using std::vector;

class DE {
 public:
 DE(std::map<std::string, std::string> mp_str, std::map<std::string, int> mp_int, std::map<std::string, double> mp_double, Problem *pro);  
 ~DE();

 virtual void run() = 0;
  double* makeIndividual();
  double getMinValID(double *pop_fitness);
  void setParentIDs(int *parent_ids, int target_id, int best_id);
  void differentialMutation(double *mutant_vector, double sf, const vector<double*> &pop, const vector<double*> &archive, int *parent_ids, int target_id, int best_id, int p_best_id);  

  string int2string(int i);
  int string2int(const string& str);
  double string2double(const string& str);
  void sortIDs(double array[], int first, int last, int index[]);
  
 protected:
  Problem *problem;

  // common parameters for DE
  string de_algorithm;
  int max_num_evaluations;
  int pop_size;
  double pop_rate;
  double scaling_factor;
  double cross_rate;
  string de_strategy;
  string de_cross;
  int max_archive_size;
  int current_archive_size;
  double archive_rate;
  double p_best_rate;
  int p_best_size;   
  //parameters for jDE
  double jde_tau_sf;
  double jde_tau_cr;
  //parameters for JADE
  double jade_learning_rate;
  //parameters for SHADE
  int shade_memory_size;

  int seed;
  string run_type;
  int run_id;
  string out_dir;

  std::mt19937 mt_seed1;
  std::mt19937 mt_seed2;

  int num_cheat_param_sampling;
  
  int dim;
  int *random_indices;
  double *mutant_vector;
  double *bsf_solution;
  double bsf_obj_val;
};

#endif
