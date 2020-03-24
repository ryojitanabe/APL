/* 
  DE.cpp
  Author: Ryoji Tanabe 2020
 */

#include "DE.h"

DE::DE(std::map<std::string, std::string> mp_str, std::map<std::string, int> mp_int, std::map<std::string, double> mp_double, Problem *pro) : problem (pro) {
    this->problem = problem;
    this->max_num_evaluations = mp_int["max_num_evaluations"];
    this->pop_size = mp_int["pop_size"];
    this->pop_rate = mp_double["pop_rate"];
    this->scaling_factor = mp_double["scaling_factor"];
    this->cross_rate = mp_double["cross_rate"];
    this->de_strategy = mp_str["de_strategy"];
    this->de_cross = mp_str["de_cross"];
    this->max_archive_size = mp_int["archive_size"];
    this->archive_rate = mp_double["archive_rate"];
    this->p_best_rate = mp_double["p_best_rate"];
    this->jde_tau_sf = mp_double["jde_tau_sf"];
    this->jde_tau_cr = mp_double["jde_tau_cr"];
    this->jade_learning_rate = mp_double["jade_learning_rate"];
    this->shade_memory_size = mp_int["shade_memory_size"];
    this->seed = mp_int["seed"];
    this->run_type = mp_str["run_type"];
    this->run_id = mp_int["run_id"];
    this->out_dir = mp_str["out_dir"];
    mt_seed1.seed(mp_int["mt_seed1"]);
    mt_seed2.seed(mp_int["mt_seed2"]);
    this->num_cheat_param_sampling = mp_int["num_cheat_param_sampling"];
    //    randUniformDouble(0,1);
    this->current_archive_size = 0;
    this->p_best_size = std::max ((int)round(pop_size * p_best_rate), 2);    
    this->dim = problem->getDimSize();

    this->random_indices = (int*)malloc(sizeof(int) * dim);
    for (int i = 0; i < dim; i++) this->random_indices[i] = i;

    this->mutant_vector = (double*)malloc(sizeof(double) * dim);
    this->bsf_solution = (double*)malloc(sizeof(double) * dim);
    this->bsf_obj_val = 1e+30;    
  }

DE::~DE(){
}

string DE::int2string(int i){
  string rt;
  std::stringstream ss;
  ss << i;
  ss >> rt;
  return rt;
}

int DE::string2int(const string& str){
  int rt;
  std::stringstream ss;
  ss << str;
  ss >> rt;
  return rt;
}

double DE::string2double(const string& str){
  double rt;
  std::stringstream ss;
  ss << str;
  ss >> rt;
  return rt;
}

 void DE::sortIDs(double array[], int first, int last, int index[]) {
   double x = array[(first + last) / 2];
   int i = first;
   int j = last;
   double temporaryValue = 0;
   int temporaryNumber = 0;

   while (true) {
     while (array[i] < x) {
       i++;
     }
     while (x < array[j]) {
       j--;
     }
     if (i >= j) {
       break;
     }

     temporaryValue = array[i];
     array[i] = array[j];
     array[j] = temporaryValue;

     temporaryNumber = index[i];
     index[i] = index[j];
     index[j] = temporaryNumber;

     i++;
     j--;
   }

   if (first < (i -1)) {
     sortIDs(array, first, i - 1, index);
   }
   if ((j + 1) < last) {
     sortIDs(array, j + 1, last, index);
   }
 }

void DE::differentialMutation(double *mutant_vector, double sf, const vector<double*> &pop, const vector<double*> &archive, int *parent_ids, int target_id, int best_id, int p_best_id) {
  int r1 = parent_ids[0];
  int r2 = parent_ids[1];
  int r3 = parent_ids[2];
  int r4 = parent_ids[3];
  int r5 = parent_ids[4];
  
  if (de_strategy == "rand_1") {
    for (int i = 0; i < dim; i++) {
      mutant_vector[i] = pop[r1][i] + sf * (pop[r2][i] - pop[r3][i]);
    }
  }
  else if (de_strategy == "rand_2") {
    for (int i = 0; i < dim; i++) {
      mutant_vector[i] = pop[r1][i] + sf * (pop[r2][i] - pop[r3][i]) + sf * (pop[r3][i] - pop[r4][i]);
    }
  }
  else if (de_strategy == "best_1") {
    for (int i = 0; i < dim; i++) {
      mutant_vector[i] = pop[best_id][i] + sf * (pop[r1][i] - pop[r2][i]);
    }
  }  
  else if (de_strategy == "best_2") {
    for (int i = 0; i < dim; i++) {
      mutant_vector[i] = pop[best_id][i] + sf * (pop[r1][i] - pop[r2][i]) + sf * (pop[r3][i] - pop[r4][i]);
    }
  }
  else if (de_strategy == "current_to_rand_1") {
    for (int i = 0; i < dim; i++) mutant_vector[i] = pop[target_id][i] + sf * (pop[r1][i] - pop[target_id][i]) + sf * (pop[r2][i] - pop[r3][i]);
  }
  else if (de_strategy == "current_to_best_1") {
    for (int i = 0; i < dim; i++) mutant_vector[i] = pop[target_id][i] + sf * (pop[best_id][i] - pop[target_id][i]) + sf * (pop[r1][i] - pop[r2][i]);
  }
  else if (de_strategy == "current_to_pbest_1") {
    if (r2 >= pop_size) {
      r2 -= pop_size;
      for (int i = 0; i < dim; i++) mutant_vector[i] = pop[target_id][i] + sf * (pop[p_best_id][i] - pop[target_id][i]) + sf * (pop[r1][i] - archive[r2][i]);
    }
    else {
      for (int i = 0; i < dim; i++) mutant_vector[i] = pop[target_id][i] + sf * (pop[p_best_id][i] - pop[target_id][i]) + sf * (pop[r1][i] - pop[r2][i]);      
    }    
  }
  else if (de_strategy == "rand_to_pbest_1") {
    if (r2 >= pop_size) {
      r2 -= pop_size;
      for (int i = 0; i < dim; i++) mutant_vector[i] = pop[r1][i] + sf * (pop[p_best_id][i] - pop[r1][i]) + sf * (pop[r2][i] - archive[r3][i]);      
    }
    else {
      for (int i = 0; i < dim; i++) mutant_vector[i] = pop[r1][i] + sf * (pop[p_best_id][i] - pop[r1][i]) + sf * (pop[r2][i] - pop[r3][i]);     
    }    
  }
  else{
    cout << "Error. " << de_strategy << " is not defined" << endl;
    exit(-1);
  }  
      
  // Repair an element out of the upper and lower bounds.
  // This method is called the "bounce-back" used in JADE
  for (int i = 0; i < dim; i++) {
    if (mutant_vector[i] < problem->getLowerBound()) mutant_vector[i] = (problem->getLowerBound() + pop[target_id][i]) / 2.0;
    if (mutant_vector[i] > problem->getUpperBound()) mutant_vector[i] = (problem->getUpperBound() + pop[target_id][i]) / 2.0;
  }  
}

void DE::setParentIDs(int *parent_ids, int target_id, int best_id) {
  std::uniform_int_distribution<int> randUniformIntPop(0, pop_size - 1);
  std::uniform_int_distribution<int> randUniformIntPopArc(0, pop_size + current_archive_size - 1);
  int r1 = -1;
  int r2 = -1;
  int r3 = -1;
  int r4 = -1;
  int r5 = -1;
  
  if (de_strategy == "rand_1" || de_strategy == "current_to_rand_1") {
    do {
      r1 = randUniformIntPop(mt_seed1);
    } while (r1 == target_id);
    do {
      r2 = randUniformIntPop(mt_seed1);
    } while ((r2 == target_id) || (r2 == r1));
    do {
      r3 = randUniformIntPop(mt_seed1);
    } while ((r3 == target_id) || (r3 == r2) || (r3 == r1));
  }
  else if (de_strategy == "rand_2") {
    do {
      r1 = randUniformIntPop(mt_seed1);
    } while (r1 == target_id);
    do {
      r2 = randUniformIntPop(mt_seed1);
    } while ((r2 == target_id) || (r2 == r1));
    do {
      r3 = randUniformIntPop(mt_seed1);
    } while ((r3 == target_id) || (r3 == r2) || (r3 == r1));
    do {
      r4 = randUniformIntPop(mt_seed1);
    } while ((r4 == target_id) || (r4 == r3) || (r4 == r2) || (r4 == r1));
    do {
      r5 = randUniformIntPop(mt_seed1);
    } while ((r5 == target_id) || (r5 == r4) || (r5 == r3) || (r5 == r2) || (r5 == r1));
  }
  else if (de_strategy == "best_1" || de_strategy == "current_to_best_1") {
    do {
      r1 = randUniformIntPop(mt_seed1);
    } while ((r1 == target_id) || (r1 == best_id));
    do {
      r2 = randUniformIntPop(mt_seed1);      
    } while ((r2 == target_id) || (r2 == best_id) || (r2 == r1));
  }
  else if (de_strategy == "best_2") {
    do {
      r1 = randUniformIntPop(mt_seed1);      
    } while ((r1 == target_id) || (r1 == best_id));
    do {
      r2 = randUniformIntPop(mt_seed1);      
    } while ((r2 == target_id) || (r2 == best_id) || (r2 == r1));
    do {
      r3 = randUniformIntPop(mt_seed1);      
    } while ((r3 == target_id) || (r3 == best_id) || (r3 == r2) ||(r3 == r1));
    do {
      r4 = randUniformIntPop(mt_seed1);      
    } while ((r4 == target_id) || (r4 == best_id) || (r4 == r3) || (r4 == r2) ||(r4 == r1));
  }
  // This implementation allows that x_r1 = x_pbest as in the original code of JADE
  else if (de_strategy == "current_to_pbest_1") {
    do {
      r1 = randUniformIntPop(mt_seed1);      
    } while ((r1 == target_id));
    do {
      r2 = randUniformIntPopArc(mt_seed1);      
    } while ((r2 == target_id) || (r2 == r1));
  }
  else if (de_strategy == "rand_to_pbest_1") {
    do {
      r1 = randUniformIntPop(mt_seed1);      
    } while ((r1 == target_id));
    do {
      r2 = randUniformIntPop(mt_seed1);      
    } while ((r2 == target_id) || (r2 == r1));
    do {
      r3 = randUniformIntPopArc(mt_seed1);      
    } while ((r3 == target_id) || (r3 == r1) || (r3 == r2));  
  }
  else{
    cout << "Error. " << de_strategy << " is not defined" << endl;
    exit(-1);
  }  
    
  parent_ids[0] = r1;
  parent_ids[1] = r2;
  parent_ids[2] = r3;
  parent_ids[3] = r4;  
  parent_ids[4] = r5;
}

double DE::getMinValID(double *pop_fitness) {
  int min_val_id = -1;
  double min_val = 1e+30;

  for (int i = 0; i < pop_size; i++) {
    if (pop_fitness[i] < min_val) {
      min_val = pop_fitness[i];
      min_val_id = i;
    }
  }

  return min_val_id;  
}

double* DE::makeIndividual() {
  std::uniform_real_distribution<double> randDoubleUniform(0.0, 1.0);
  double* x = (double*)malloc(sizeof(double) * dim);

  for (int i = 0; i < dim; i++) {
    x[i] = ((problem->getUpperBound() - problem->getLowerBound()) * randDoubleUniform(mt_seed1)) + problem->getLowerBound();
  }

  return x;
}
