/* 
  Problem.cpp
  Author: Ryoji Tanabe 2020
 */

#include "Problem.h"
#include "bbobStructures.h"

#include<iostream>
#include <sstream>
#include<fstream>

Problem::Problem(string test_func, int dim, string variable_order, string func_type, int bbob_func_id) {
    this->dim = dim;
    this->variable_order = variable_order;
    this->func_type = func_type;
    
    if (func_type == "bbob") this->test_func = "bbob-f" + int2string(bbob_func_id);
    else this->test_func = test_func;

    //   initializeRandomFuncitons();
//   eval_vec = (variable*)malloc(sizeof(variable) * problem_size);
//   shifted_eval_vec = (variable*)malloc(sizeof(variable) * problem_size);
//   r_temp_vec = (variable*)malloc(sizeof(variable) * 2);
//   temp_vec_for_trans = (variable*)malloc(sizeof(variable) * 2);    
  }

Problem::~Problem() {}

string Problem::int2string(int i) {
  string rt;
  std::stringstream ss;
  ss << i;
  ss >> rt;
  return rt;
}

int Problem::string2int(const string& str) {
  int rt;
  std::stringstream ss;
  ss << str;
  ss >> rt;
  return rt;
}

double Problem::string2double(const string& str) {
  double rt;
  std::stringstream ss;
  ss << str;
  ss >> rt;
  return rt;
}

double Problem::errorValBSF() {
  if (fabs(bsf_obj_val - opt_obj_val) < val2rearch) return 0;
  else return fabs(bsf_obj_val - opt_obj_val);
}

bool Problem::updateBSFsolution(double f, double *x) {
  if (f < bsf_obj_val) {
    bsf_obj_val = f;
    for (int i = 0; i < dim; i++) bsf_solution[i] = x[i];    
    return true;
  }
  else {
    return false;
  }  
}

void Problem::setInstanceData() {
  if (func_type == "bbob") {
    val2rearch = 1e-8;
    opt_obj_val = fgeneric_ftarget();
    bsf_obj_val = 1e+30;
    bsf_solution = new double[dim];

    upper_bound = 5;
    lower_bound = -5;    
  }
  else {
    val2rearch = 1e-8;
    opt_obj_val = 0;
    bsf_obj_val = 1e+30;
    bsf_solution = new double[dim];
    shifted_opt_pos = new double[dim];
    random_permutation = new int[dim];
    tmp_vec = new double[dim];
    translated_x = new double[dim];
  
    string in_file;
    string temp_str;
  
    in_file = "func_data/" + test_func + "_shifted_data.dat";
    std::ifstream ifs_s(in_file.c_str());

    if(!ifs_s) {
      cout <<"Error. " <<  in_file << " does not exist" << endl;
      exit(-1);
    }

    for (int i = 0; i < dim; i++) {  
      ifs_s >> temp_str;
      shifted_opt_pos[i] = string2double(temp_str);
      shifted_opt_pos[i] = 0;
    }

    ifs_s.close();

    in_file = "func_data/rand_perm/D" + int2string(dim) + ".dat";
    std::ifstream ifs_rp(in_file.c_str());

    if (!ifs_rp) {
      cout <<"Error. " <<  in_file << " does not exist" << endl;
      exit(-1);
    }

    if (variable_order == "alphanumeric") {
      for (int i = 0; i < dim; i++) random_permutation[i] = i;  
    }
    else if(variable_order == "random") {
      for (int i = 0; i < dim; i++) {  
	ifs_rp >> temp_str;
	random_permutation[i] = string2int(temp_str);
      }
    }
    else {
      cout <<"Error. " <<  variable_order << " is not defined" << endl;
      exit(-1);
    }

    ifs_rp.close();

    if (test_func== "brock-rotated-ellipsoid" || test_func == "brock-rotated-ackley") {
      d2_rot_matrix = new double**[dim];
      for (int i = 0; i < dim; i++) {
	d2_rot_matrix[i] = new double*[2];
	for (int j = 0; j < 2; j++) d2_rot_matrix[i][j] = new double[2];
      }

      for (int i= 0; i < dim; i++) {
	in_file = "func_data/rotation/rot_matrix/D2_" + int2string(i) + ".dat";
	std::ifstream ifs_rmat(in_file.c_str());

	if (!ifs_rp) {
	  cout <<"Error. " <<  in_file << " does not exist" << endl;
	  exit(-1);
	}

	for (int j = 0; j < 2; j++) {  
	  for (int k = 0; k < 2; k++) {  
	    ifs_rmat >> temp_str;
	    d2_rot_matrix[i][j][k] = string2double(temp_str);     
	  }
	}

	// for (int j = 0; j < 2; j++) {  
	// 	for (int k = 0; k < 2; k++) {  
	// 	  cout << d2_rot_matrix[i][j][k] <<" ";
	// 	}
	// 	cout << endl;
	// }

	ifs_rmat.close();
      }
    }

    if (test_func== "fully-rotated-ellipsoid" || test_func== "fully-rotated-step-ellipsoid") {
      rot_matrix = (double **)malloc(sizeof(double) * dim);
      for (int i = 0; i < dim; i++) rot_matrix[i] = (double *)malloc(sizeof(double) * dim);
    
      in_file = "func_data/rotation/rot_matrix/D" + int2string(dim) + ".dat";
      std::ifstream ifs_m(in_file.c_str());

      if (!ifs_m) {
	cout <<"Error. " <<  in_file << " does not exist" << endl;
	exit(-1);
      }

      for (int i = 0; i < dim; i++) {  
	for (int j = 0; j < dim; j++) {  
	  ifs_m >> temp_str;
	  rot_matrix[i][j] = string2double(temp_str);
	}
      }

      ifs_m.close();
    }

  
    // cout <<"---------------------------" << endl;
    // for (int i = 0; i < dim; i++) {  
    // 	for (int j = 0; j < dim; j++) {  
    // 	  cout << rot_matrix[i][j] << " ";
    // 	}
    // 	cout << endl;
    // }
    // cout <<"---------------------------" << endl;

    if (test_func== "sphere") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else if (test_func== "schwefels_222") {
      upper_bound = 10;
      lower_bound = -10;
    }
    else if (test_func== "schwefels_12") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else if (test_func== "schwefels_221") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else if (test_func== "rosenbrock") {
      upper_bound = 30;
      lower_bound = -30;
    }
    else if (test_func== "step") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else if (test_func== "quartic_noize") {
      upper_bound = 1.28;
      lower_bound = -1.28;
      val2rearch = pow(10.0, -2);
    }
    else if (test_func== "schwefels_226") {
      upper_bound = 500;
      lower_bound = -500;
    }
    else if (test_func== "rastrigin") {
      upper_bound = 5.12;
      lower_bound = -5.12;
    }
    else if (test_func== "ackley" || test_func == "t2_schwefels_221_ackley")  {
      upper_bound = 32.0;
      lower_bound = -32.0;
    }
    else if (test_func== "griewank") {
      upper_bound = 600.0;
      lower_bound = -600.0;
    }
    else if (test_func== "t-griewank") {
      upper_bound = 600.0;
      lower_bound = -600.0;
    }
    else if (test_func== "penalized1") {
      upper_bound = 50.0;
      lower_bound = -50.0;
    }
    else if (test_func== "penalized2") {
      upper_bound = 50.0;
      lower_bound = -50.0;
    }
    else if (test_func== "schafferF7") {
      upper_bound = 10;
      lower_bound = -10;
    }
    else if (test_func== "weierstrass") {
      upper_bound = 5.0;
      lower_bound = -5.0;
    }
    else if (test_func== "katsuura") {
      upper_bound = 0.5;
      lower_bound = -0.5;
    }
    else if (test_func== "t2_schwefels_221" || test_func== "t3_schwefels_221" || test_func== "t4_schwefels_221") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else if (test_func== "t_hf_rastrigin") {
      // upper_bound = 5.12;
      // lower_bound = -5.12;
      upper_bound = 32;
      lower_bound = -32;
    }
    else if (test_func== "exp_gri_ros") {
      upper_bound = 30;
      lower_bound = -30;
    }
    else if (test_func== "brock-rotated-ackley") {
      upper_bound = 32;
      lower_bound = -32;
    }
    else if (test_func== "ellipsoid" || test_func== "brock-rotated-ellipsoid" || test_func== "fully-rotated-ellipsoid" || test_func== "fully-rotated-step-ellipsoid" || test_func== "step-ellipsoid") {
      upper_bound = 5;
      lower_bound = -5;
    }
    else if (test_func== "bohachevsky") {
      upper_bound = 100;
      lower_bound = -100;
    }
    else {
      cout << "Error. " <<test_func << " is not defined" << endl;
      exit(-1);
    }        
  }  
}

double Problem::eval(double *x) {
  double f = 0;
  //  for (int i = 0; i < dim; i++) eval_vec[i] = individual[rand_order_vec[i]] - rand_shifted_vec[rand_order_vec[i]];
  // for (int i = 0; i < dim; i++) eval_vec[i] = individual[rand_order_vec[i]] - rand_shifted_vec[i];

  if (func_type == "bbob") {
    f = fgeneric_evaluate(x);
  }
  else { 
    for (int i = 0; i < dim; i++) tmp_vec[i] = x[i] - shifted_opt_pos[i];
    for (int i = 0; i < dim; i++) translated_x[i] = tmp_vec[random_permutation[i]];  

    if (test_func== "sphere") f = sphere(translated_x);
    else if (test_func== "rosenbrock") f = rosenbrock(translated_x);
    else if (test_func== "rastrigin") f = rastrigin(translated_x);
    else {
      cout <<"Error. " << test_func << "is not defined" << endl;
      exit(-1);
    }

    //  if (fabs(f - opt_obj_val) < val2rearch) f = 0;
  }
  
  return f;
}

double Problem::sphere(double *x) {
  double f = 0;
  for (int i = 0; i < dim; i++) f += x[i] * x[i];
  return f;
}

double Problem::rosenbrock(double *x) {
  double f = 0;
  for (int i = 0; i < dim - 1; i++) f += (100 * pow((x[i] * x[i] - x[i + 1]), 2.0)) + pow((x[i] - 1.0), 2.0);  
  return f;
}

double Problem::rastrigin(double *x) {
  double f = 10 * dim;
  for (int i = 0; i < dim; i++) f += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);    
  return f;
}

