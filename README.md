# Analyzing Adaptive Parameter Landscapes in Parameter Adaptation Methods for Differential Evolution 

This repository provides the C++ code used in the following paper:

> Ryoji Tanabe, **Analyzing Adaptive Parameter Landscapes in Parameter Adaptation Methods for Differential Evolution**, ACM GECCO2020, [pdf](https://ryojitanabe.github.io/pdf/t-gecco2020.pdf)

This code is based on the *old* version of the COCO framework implemented in C (https://coco.gforge.inria.fr/):

> Nikolaus Hansen, Anne Auger, Olaf Mersmann, Tea Tusar, Dimo Brockhoff, **COCO: A Platform for Comparing Continuous Optimizers in a Black-Box Setting**, CoRR abs/1603.08785 (2016), [link](https://arxiv.org/abs/1603.08785)

# Requirements

I compiled this code with gcc version 7.5.0 on the Ubuntu 18.04 LTS. Python (2 or 3) is needed to run "src/createfolders.py". Ruby (>=2.5) is also needed to run "src/run_de.rb", which is an optional script to perform an experiment automatically.

# Usage

## Compile

```
$ git clone https://github.com/ryojitanabe/APL.git
$ cd APL
$ make
```

## Simple example

```
$ python src/createfolders.py tmp 
$ ./de -coco_func_num 1 -prosize 20 -nfe 200000 -alg shade -pop 100 -strategy current_to_pbest_1 -p_rate 0.05 -cross bin -pop_rate -1 -arc_rate 1 -memory_size 10 -num_cheat_param_sampling 50 -out_dir tmp -run_type bbob-apl-analysis -in_file dummy.dat
```

The first Python command creates a directory "tmp" for the procedure of the COCO framework. The second command runs a DE algorithm on the 20-dimensional 3rd function instance of f1 (the Sphere function). More precisely, the parameters of the DE are as follows:

- Maximum number of function evaluations (nfe): 200000 (=10,000 * 20)
- Parameter adaptation method (alg): SHADE
- Population size (pop): 100
- Mutation strategy (strategy): current-to-pbest/1
- Pbest rate in current-to-pbest/1 (p\_rate): 0.05
- Crossover method (cross): binomial crossover
- Population rate (pop\_rate): -1
    - This sets the population size to "dim * pop\_rate".
	- "pop\rate_=-1" means that the population size is determined by the argument "pop". 
- Archive rate (arc\_rate): 1
    - The actual archive size is pop * arc\_rate. Thus, in this example, 100 * 1 = 100.
- Memory size in SHADE (memory\_size): 10
- Memory size in SHADE (num\_cheat\_param\_sampling): 50
    - This generates 50 * 50 parameter pairs of the scale factor F and the crossover rate C.
- Output directory (out_dir):  tmp
- Run type (run\_type): bbob-apl-analysis
    - This runs DE to obtain adaptive parameter landscapes.
	- When running DE for benchmarking, it should be set to "bbob".
- Problem instance (in\_file): dummy.dat
	- This determines which function instance should be used for the experiment.

This command outputs many dat files in the "tmp" directory such as "parameters\_bbob-f1\_d20\_9th\_ranked\_ind\_11000th\_evals\_3th\_run.dat". This file contains 50 * 50 parameter pairs (F and C) and their corresponding fitness improvement values (the G1 value) of the 9th *ranked* individual in the population at 1,000 function evaluations on the 20-dimensional 3rd function instance of f1.
 
Similarly, the following is the command to run a DE with the parameter adaptation method in jDE:

```
$ python src/createfolders.py tmp 
$ ./de -coco_func_num 1 -prosize 20 -nfe 200000 -alg jde -pop 100 -strategy current_to_pbest_1 -p_rate 0.05 -cross bin -pop_rate -1 -arc_rate 1 -tau_sf 0.1 -tau_cr 0.1 -num_cheat_param_sampling 50 -out_dir tmp -run_type bbob-apl-analysis -in_file dummy.dat
```

The following is also the command to run a DE with the parameter adaptation method in JADE:

```
$ python src/createfolders.py tmp 
$ ./de -coco_func_num 1 -prosize 20 -nfe 200000 -alg jade -pop 100 -strategy current_to_pbest_1 -p_rate 0.05 -cross bin -pop_rate -1 -arc_rate 1 -jade_learning_rate 0.1 -num_cheat_param_sampling 50 -out_dir tmp -run_type bbob-apl-analysis -in_file dummy.dat
```

## Reproduce results presented in the GECCO paper

The Ruby script "run_de.rb" automatically performs all experiments to reproduce all results presented in the GECCO paper:

```
$ ruby run_de.rb
```

This command outputs results in the "results" directory. The Ruby script refers to a file in "median_run_ids", which provides an index number of a single run with a median performance among 15 runs for each function.
