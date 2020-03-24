/* 
  ClassicalDE.cpp
  Author: Ryoji Tanabe 2020
 */

#include "ClassicalDE.h"

void classicalDE::run() {
  cout <<  std::scientific << std::setprecision(8);
  std::uniform_real_distribution<double> randUniformDouble(0.0, 1.0);
  std::uniform_int_distribution<int> randUniformIntPop(0, pop_size - 1);
  std::uniform_int_distribution<int> randUniformIntDim(0, dim - 1);
  std::uniform_int_distribution<int> randUniformIntPbest(0, p_best_size - 1);
  std::uniform_int_distribution<int> randUniformIntArchive(0, max_archive_size - 1);

  int nfes = 0;
  double error_val = 1e+30;

  std::string profile_out_file = out_dir + "/" + int2string(run_id) + "th_run.dat";
  std::ofstream profile_ofs(profile_out_file);
  //  profile_ofs << std::setprecision(5);
  profile_ofs <<  std::scientific << std::setprecision(5);  

  //initialize population
  std::vector <double*> pop;
  double *pop_fitness = (double*)malloc(sizeof(double) * pop_size);
  std::vector <double*> children;
  double *children_fitness = (double*)malloc(sizeof(double) * pop_size);
  std::vector <double*> archive;
  
  for (int i = 0; i < pop_size; i++) {
    pop.push_back((double*)malloc(sizeof(double) * dim));
    children.push_back((double*)malloc(sizeof(double) * dim));
  }

  for (int i = 0; i < max_archive_size; i++) {
    archive.push_back((double*)malloc(sizeof(double) * dim));
  }
  
  for (int i = 0; i < pop_size; i++) {
    pop[i] = makeIndividual();
    //for (int j = 0; j < dim; j++) pop[i][j] = ((max_region - min_region) * randUniformDouble(mt_seed1)) + min_region;
  }

  for (int i = 0; i < pop_size; i++) {
    pop_fitness[i] = problem->eval(pop[i]);
    nfes++;
    if (problem->updateBSFsolution(pop_fitness[i], pop[i]) == true) {
	error_val = problem->errorValBSF();
	//	cout << nfes << " " << error_val << endl;
	profile_ofs << nfes << " " << error_val << endl;
	//	if (error_val <= 0) exit(0);
	if (error_val <= 0) break;
    }
    
    if (nfes >= max_num_evaluations) break;
  }
  
  int *parent_ids = (int*)malloc(sizeof(int) * 5);
  int best_id;
  int p_best_id;
  int rnd_id;
  double *mutant_vector = (double*)malloc(sizeof(double) * dim);
  int cross_rand_pos;
  double *cross_rnd_values = (double*)malloc(sizeof(double) * dim);

  int *sorted_ids = (int*)malloc(sizeof(int) * pop_size);
  double *tmp_vals = (double*)malloc(sizeof(double) * pop_size);

  // for the parameter landscape analysis    
  double *virtual_child = (double*)malloc(sizeof(double) * dim);  
  double virtual_child_fitness;
  double fitness_improvement;
  int each_row_column_size = num_cheat_param_sampling;
  double *gird_sf_vals = (double*)malloc(sizeof(double) * each_row_column_size);  
  double *gird_cr_vals = (double*)malloc(sizeof(double) * each_row_column_size);  
  double grid_val = 0.0;
  double tmp_sf_val;
  double tmp_cr_val;
  int *individual_ranks = (int*)malloc(sizeof(int) * pop_size);  
  int count_iter = 1;
  
  for (int j = 0; j < each_row_column_size; j++) {
    gird_sf_vals[j] = grid_val;
    gird_cr_vals[j] = grid_val;
    grid_val += 1.0 / (each_row_column_size - 1);
  }
  
  while (nfes < max_num_evaluations && error_val > 0) {
    best_id = getMinValID(pop_fitness);
    
    for (int i = 0; i < pop_size; i++) sorted_ids[i] = i;
    for (int i = 0; i < pop_size; i++) tmp_vals[i] = pop_fitness[i];
    sortIDs(&tmp_vals[0], 0, pop_size - 1, sorted_ids);

    // for the parameter landscape analysis     
    for (int i = 0; i < pop_size; i++) individual_ranks[sorted_ids[i]] = i;
    
    for (int i = 0; i < pop_size; i++) {
      p_best_id = sorted_ids[randUniformIntPbest(mt_seed1)];
      
      setParentIDs(parent_ids, i, best_id);      
      differentialMutation(mutant_vector, scaling_factor, pop, archive, parent_ids, i, best_id, p_best_id);
      
      //      crossover(children[i], mutant_vector, pop[i]);
      //crossover
      cross_rand_pos = randUniformIntDim(mt_seed1);
      for (int j = 0; j < dim; j++) cross_rnd_values[j] = randUniformDouble(mt_seed1);
	
      for (int j = 0; j < dim; j++) {
  	if (cross_rnd_values[j] <= cross_rate || j == cross_rand_pos) children[i][j] = mutant_vector[j];
  	else children[i][j] = pop[i][j];
      }

      if (run_type == "bbob-apl-analysis") {
	// for the parameter landscape analysis
	if (count_iter % 10 == 0 || count_iter == 1) {
	  //	  cout << count_iter << endl;
	  std::string la_out_file = out_dir + "/parameters_" + problem->getTestFuncName() + "_d" + int2string(dim) + "_" + int2string(individual_ranks[i]) + "th_ranked_ind_" + int2string(nfes) + "th_evals_" + int2string(run_id) + "th_run.dat";
	  std::ofstream ofs(la_out_file);
	  ofs <<  std::setprecision(5);

	  for (int j = 0; j < each_row_column_size; j++) {
	    for (int k = 0; k < each_row_column_size; k++) {
	      tmp_sf_val = gird_sf_vals[j];
	      tmp_cr_val = gird_cr_vals[k];
	  
	      // mutation
	      differentialMutation(mutant_vector, tmp_sf_val, pop, archive, parent_ids, i, best_id, p_best_id);
	      // crossover	
	      for (int l = 0; l < dim; l++) {
		if (cross_rnd_values[l] <= tmp_cr_val || l == cross_rand_pos) virtual_child[l] = mutant_vector[l];
		else virtual_child[l] = pop[i][l];
	      }

	      virtual_child_fitness = problem->eval(virtual_child);
	      fitness_improvement = pop_fitness[i] - virtual_child_fitness;
	      fitness_improvement = std::max(fitness_improvement, 0.0);
	      ofs << fitness_improvement << endl;
	    }
	  }
	  
	  ofs.close();
	}
	// for the parameter landscape analysis
		
      }      
    }

    for (int i = 0; i < pop_size; i++) {
      children_fitness[i] = problem->eval(children[i]);
      nfes++;
      if (problem->updateBSFsolution(children_fitness[i], children[i]) == true) {
	error_val = problem->errorValBSF();
	//	cout << nfes << " " << error_val << endl;
	profile_ofs << nfes << " " << error_val << endl;
	//	if (error_val <= 0) exit(0);
	if (error_val <= 0) break;
      }

      if (nfes >= max_num_evaluations) break;
    }

    for (int i = 0; i < pop_size; i++) {
      if (children_fitness[i] <= pop_fitness[i]) {
	if (current_archive_size < max_archive_size) {
	  for (int j = 0; j < dim; j++) archive[current_archive_size][j] = pop[i][j];
	      current_archive_size++;
	}
	else {	  
	  rnd_id = randUniformIntArchive(mt_seed1);
	  for (int j = 0; j < dim; j++) archive[rnd_id][j] = pop[i][j];
	}

	pop_fitness[i] = children_fitness[i];
	for (int j = 0; j < dim; j ++)  pop[i][j] = children[i][j];
      }
    }

    count_iter++;    
  }

  profile_ofs.close();  
}
