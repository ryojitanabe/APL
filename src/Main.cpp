/* 
  Main.cpp

  Author: Ryoji Tanabe 2020
 */

#include "ClassicalDE.h"
#include "JDE.h"
#include "JADE.h"
#include "SHADE.h"

#include "bbobStructures.h" /* Include all declarations for BBOB calls */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <map>

// parameters for test functions
int g_coco_func_num;
string g_test_func;
int g_dim;
string g_variable_order;
// common parameters for DE
string g_de_algorithm;
int g_max_num_evaluations;
int g_pop_size;
double g_pop_rate = -1;
double g_scaling_factor;
double g_cross_rate;
string g_de_strategy;
string g_de_cross;
int g_archive_size;
double g_archive_rate = -1;
double g_p_best_rate;
int g_seed;
string g_run_type;
int g_run_id;
string g_out_dir;
//parameters for jDE
double g_jde_tau_sf;
double g_jde_tau_cr;
//parameters for JADE
double g_jade_learning_rate;
//parameters for SHADE
int g_shade_memory_size;

// parameters for analyzing adaptive parameter landscapes
int g_mt_seed1;
int g_mt_seed2;
int g_num_cheat_param_sampling;
string g_in_file;

int parseCommandLine(int argc, char **argv) {
  int index = 1;
  char *flag, *value;

  while (index < argc) { // read (flab,value) pairs from command line arguments
    flag = argv[index];
    value = argv[index+1];
    //    cout << flag << " " << value << endl;

    if (!strcmp(flag, "-test_func")) {
      g_test_func = value;
    }
    else if (!strcmp(flag, "-coco_func_num")) {
      g_coco_func_num = atoi(value);
    }
    else if (!strcmp(flag, "-variable_order")) {
      g_variable_order = value;
    }
    else if (!strcmp(flag, "-prosize")) {
      g_dim = atoi(value);
    }
    else if (!strcmp(flag, "-nfe")) {
      g_max_num_evaluations = atoi(value);
    }
    else if (!strcmp(flag, "-alg")) {
      g_de_algorithm = value;
    }
    else if (!strcmp(flag,"-pop_rate")) {
      g_pop_rate = atof(value);
    }
    else if (!strcmp(flag,"-arc_rate")) {
      g_archive_rate = atof(value);
    }
    else if (!strcmp(flag,"-pop")) {
      g_pop_size = atoi(value);
    }
    else if (!strcmp(flag,"-strategy")) {
      g_de_strategy = value;
    }
    else if (!strcmp(flag,"-cross")) {
      g_de_cross = value;
    }
    else if (!strcmp(flag,"-sf")) {
      g_scaling_factor = atof(value);
    }
    else if (!strcmp(flag,"-cr")) {
      g_cross_rate = atof(value);
    }
    else if (!strcmp(flag, "-arc_size")) {
      g_archive_size = atoi(value);
    }
    else if (!strcmp(flag, "-memory_size")) {
      g_shade_memory_size = atoi(value);
    }
    else if (!strcmp(flag, "-jade_learning_rate")) {
      g_jade_learning_rate = atof(value);
    }
    else if (!strcmp(flag, "-p_rate")) {
      g_p_best_rate = atof(value);
    }
    else if (!strcmp(flag, "-tau_sf")) {
      g_jde_tau_sf = atof(value);
    }
    else if (!strcmp(flag, "-tau_cr")) {
      g_jde_tau_cr = atof(value);
    }
    else if (!strcmp(flag, "-run_type")) {
      g_run_type = value;
    }
    else if (!strcmp(flag, "-run_id")) {
      g_run_id = atoi(value);
    }
    else if (!strcmp(flag, "-out_dir")) {
      g_out_dir = value;
    }
    // else if (!strcmp(flag, "-mt_seed1")) {
    //   g_mt_seed1 = atoi(value);
    // }
    // else if (!strcmp(flag, "-mt_seed2")) {
    //   g_mt_seed2 = atoi(value);
    // }
    else if (!strcmp(flag, "-num_cheat_param_sampling")) {
      g_num_cheat_param_sampling = atoi(value);
    }
    else if (!strcmp(flag, "-in_file")) {
      g_in_file = value;
    }
    else if (!strcmp(flag, "-seed")) {
      g_seed = atoi(value);
      srand(g_seed);
    }
    else {
      cout << "Error. Unkown argument " << flag << endl;
      exit(-1);
    }

    index += 2;
  }

  return index;
}

int main(int argc, char **argv) {
  srand(0);
  parseCommandLine(argc, argv);
    
  unsigned int instances[15] = {1, 2, 3, 4, 5, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
  unsigned int idx_instances, seed;
  int independent_restarts;
  double maxfunevals, minfunevals;

  clock_t t0 = clock();
  time_t Tval;
  /**************************************************
   *          BBOB Mandatory initialization         *
   *************************************************/
  /* retrieve all default parameters of BBOB calls  */
  ParamStruct params = fgeneric_getDefaultPARAMS();

  /* modify the following parameters, choosing a different setting
   * for each new experiment */
  strcpy(params.dataPath, g_out_dir.c_str());  /* different folder for each experiment! */
  /* please beforehand run from the command-line 'python createfolders.py PUT_MY_BBOB_DATA_PATH'
   * to create the necessary folder structure to run an experiment. */
  strcpy(params.algName, "g_de_algorithm");
  strcpy(params.comments, "PUT MORE DETAILED INFORMATION, PARAMETER SETTINGS ETC");

  seed = time(NULL);
  srand(seed); /* used by MY_OPTIMIZER */
  printf("random seed set to %d\n", seed);

  /* To make the noise deterministic. */
  /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */

  /* now the main loop */

  // This performs DE only on a single problem instance for each function and each dimensionality
  if (g_run_type == "bbob-apl-analysis" || g_run_type == "bbob-adaptation") {    
    std::ifstream ifs(g_in_file);
    if (!ifs) {
      cout << "Error. " + g_in_file + " does not exist" << endl;
      cout << "This program needs the file that specifies a single run with a median best-so-far error value among 15 runs" << endl;
      exit(-1);     
    }

    std::string str;
    ifs >> str;    
    idx_instances = std::stoi(str);
    
    /* set DIM, funcId, instanceId to initialize BBOB fgeneric */
    params.DIM = g_dim;
    params.funcId = g_coco_func_num;
    params.instanceId = instances[idx_instances];
    /* call the BBOB initialization */
    fgeneric_initialize(params);

    /* 5. * dim should be fine to just check everything */
    minfunevals = 0;//g_problem_size + 2;  /* PUT MINIMAL USEFUL NUMBER OF FEVALS */
    independent_restarts = -1;

    g_run_id = idx_instances;
    g_mt_seed1 = g_run_id;
    g_mt_seed2 = g_run_id + 100;
    
    if (g_dim > 3 && g_pop_rate > 0) g_pop_size = (int)round(g_pop_rate * g_dim);
    else if (g_pop_rate > 0) g_pop_size = 20;
    if (g_archive_rate > 0) g_archive_size = (int)round(g_pop_size * g_archive_rate);

    std::map<std::string, std::string> mp_str;
    std::map<std::string, int> mp_int;
    std::map<std::string, double> mp_double;
  
    mp_int["max_num_evaluations"] = g_max_num_evaluations;
    mp_int["pop_size"] = g_pop_size;
    mp_double["pop_rate"] = g_pop_rate;
    mp_double["scaling_factor"] = g_scaling_factor;
    mp_double["cross_rate"] = g_cross_rate;
    mp_str["de_strategy"] = g_de_strategy;
    mp_str["de_cross"] = g_de_cross;
    mp_int["archive_size"] = g_archive_size;
    mp_double["archive_rate"] = g_archive_rate;
    mp_double["p_best_rate"] = g_p_best_rate;
    mp_double["jde_tau_sf"] = g_jde_tau_sf;
    mp_double["jde_tau_cr"] = g_jde_tau_cr;
    mp_double["jade_learning_rate"] = g_jade_learning_rate;
    mp_int["shade_memory_size"] = g_shade_memory_size;
    mp_int["seed"] = g_seed;
    mp_str["run_type"] = g_run_type;
    mp_int["run_id"] = g_run_id;
    mp_str["out_dir"] = g_out_dir;
    mp_int["mt_seed1"] = g_mt_seed1;
    mp_int["mt_seed2"] = g_mt_seed2;
    mp_int["num_cheat_param_sampling"] = g_num_cheat_param_sampling;
    
    Problem problem(g_test_func, g_dim, g_variable_order, "bbob", g_coco_func_num);
    problem.setInstanceData();
    
    DE *alg;  
    if (g_de_algorithm == "classical-de") alg = new classicalDE(mp_str, mp_int, mp_double, &problem);
    else if (g_de_algorithm == "jde") alg = new JDE(mp_str, mp_int, mp_double, &problem);
    else if (g_de_algorithm == "jade") alg = new JADE(mp_str, mp_int, mp_double, &problem);
    else if (g_de_algorithm == "shade") alg = new SHADE(mp_str, mp_int, mp_double, &problem);
    else {
      cout << "Error. " << g_de_algorithm << " is not defined" << endl;
      exit(-1);
    }    
    
    alg->run();
    delete alg;
      
    printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", g_coco_func_num, g_dim,
	   instances[idx_instances], fgeneric_evaluations(), independent_restarts);
    printf(" fbest-ftarget=%.4e, elapsed time [h]: %.2f\n", 
	   fgeneric_best() - fgeneric_ftarget(), (double)(clock()-t0)/CLOCKS_PER_SEC/60./60.);
    /* call the BBOB closing function to wrap things up neatly */
    fgeneric_finalize();    
  }
  else {
    // This performs DE on all the 15 problem instances for each function and each dimensionality.
    // Thus, this is the traditional procedure of BBOB.
    
    /* Function indices are from 1 to 24 (noiseless) or from 101 to 130 (noisy) */
    /* for the noisy functions exchange the for loop with */
    /* for (ifun = 101; ifun <= 130; ifun++) */
    for (idx_instances = 0; idx_instances < 15; idx_instances++) {
      /* set DIM, funcId, instanceId to initialize BBOB fgeneric */
      params.DIM = g_dim;
      params.funcId = g_coco_func_num;
      params.instanceId = instances[idx_instances];
      /* call the BBOB initialization */
      fgeneric_initialize(params);

      /* 5. * dim should be fine to just check everything */
      minfunevals = 0;//g_problem_size + 2;  /* PUT MINIMAL USEFUL NUMBER OF FEVALS */
      independent_restarts = -1;

      g_run_id = idx_instances;
      g_mt_seed1 = g_run_id;
      g_mt_seed2 = g_run_id + 100;
    
      if (g_dim > 3 && g_pop_rate > 0) g_pop_size = (int)round(g_pop_rate * g_dim);
      else if (g_pop_rate > 0) g_pop_size = 20;
      if (g_archive_rate > 0) g_archive_size = (int)round(g_pop_size * g_archive_rate);

      std::map<std::string, std::string> mp_str;
      std::map<std::string, int> mp_int;
      std::map<std::string, double> mp_double;
  
      mp_int["max_num_evaluations"] = g_max_num_evaluations;
      mp_int["pop_size"] = g_pop_size;
      mp_double["pop_rate"] = g_pop_rate;
      mp_double["scaling_factor"] = g_scaling_factor;
      mp_double["cross_rate"] = g_cross_rate;
      mp_str["de_strategy"] = g_de_strategy;
      mp_str["de_cross"] = g_de_cross;
      mp_int["archive_size"] = g_archive_size;
      mp_double["archive_rate"] = g_archive_rate;
      mp_double["p_best_rate"] = g_p_best_rate;
      mp_double["jde_tau_sf"] = g_jde_tau_sf;
      mp_double["jde_tau_cr"] = g_jde_tau_cr;
      mp_double["jade_learning_rate"] = g_jade_learning_rate;
      mp_int["shade_memory_size"] = g_shade_memory_size;
      mp_int["seed"] = g_seed;
      mp_str["run_type"] = g_run_type;
      mp_int["run_id"] = g_run_id;
      mp_str["out_dir"] = g_out_dir;
      mp_int["mt_seed1"] = g_mt_seed1;
      mp_int["mt_seed2"] = g_mt_seed2;
      mp_int["num_cheat_param_sampling"] = g_num_cheat_param_sampling;
    
      Problem problem(g_test_func, g_dim, g_variable_order, "bbob", g_coco_func_num);
      problem.setInstanceData();
    
      DE *alg;  
      if (g_de_algorithm == "classical-de") alg = new classicalDE(mp_str, mp_int, mp_double, &problem);
      else if (g_de_algorithm == "jde") alg = new JDE(mp_str, mp_int, mp_double, &problem);
      else if (g_de_algorithm == "jade") alg = new JADE(mp_str, mp_int, mp_double, &problem);
      else if (g_de_algorithm == "shade") alg = new SHADE(mp_str, mp_int, mp_double, &problem);
      else {
	cout << "Error. " << g_de_algorithm << " is not defined" << endl;
	exit(-1);
      }    
    
      alg->run();
      delete alg;
      
      printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", g_coco_func_num, g_dim,
	     instances[idx_instances], fgeneric_evaluations(), independent_restarts);
      printf(" fbest-ftarget=%.4e, elapsed time [h]: %.2f\n", 
	     fgeneric_best() - fgeneric_ftarget(), (double)(clock()-t0)/CLOCKS_PER_SEC/60./60.);
      /* call the BBOB closing function to wrap things up neatly */
      fgeneric_finalize();
    }
  }
  
  Tval = time(NULL);

  return 0;
}
