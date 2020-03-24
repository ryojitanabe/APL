# -*- coding: utf-8 -*-
include Math

pop_rate = -1
pop_size = 100
scaling_factor = 0.5
cross_rate =  0.9
jade_learning_rate = 0.1
tau_sf = 0.1
tau_cr = 0.1
memory_size = 10
p_best_rate = 0.05
arc_rate = 1
num_cheat_param_sampling = 50
de_strategy = "current_to_pbest_1"
de_cross = "bin"

run_type = "bbob-apl-analysis"
["jde", "jade", "shade"].each{|de_algorithm|
  out_dir = "results/pal_#{de_algorithm}_n#{pop_size}_#{de_strategy}_#{de_cross}" 

  if !File.directory?(out_dir)
    Dir.mkdir(out_dir)
  end

  [2, 3, 5, 10, 20, 40].each{|problem_size|
    max_num_evaluations = problem_size * (10**4)  
    (1..24).each{|coco_func_num|          
      instance_out_dir = "#{out_dir}/f#{coco_func_num}_d#{problem_size}"
      puts instance_out_dir
      
      if !File.directory?(instance_out_dir)
        Dir.mkdir(instance_out_dir)
      end

      system("python src/createfolders.py #{out_dir}")
      in_file = "median_run_ids/#{de_algorithm}_n#{pop_size}_#{de_strategy}_#{de_cross}/f#{coco_func_num}_d#{problem_size}/median_run_id.dat"
      system("./de -coco_func_num #{coco_func_num} -prosize #{problem_size} -nfe #{max_num_evaluations} -alg #{de_algorithm} -pop #{pop_size} -strategy #{de_strategy} -cross #{de_cross} -sf #{scaling_factor} -cr #{cross_rate} -pop_rate #{pop_rate} -arc_rate #{arc_rate}  -memory_size #{memory_size} -jade_learning_rate #{jade_learning_rate} -p_rate #{p_best_rate}  -tau_sf #{tau_sf} -tau_cr #{tau_cr} -num_cheat_param_sampling #{num_cheat_param_sampling} -out_dir #{instance_out_dir} -run_type #{run_type} -in_file #{in_file}")
    }
  }
}


# # For benchmarking DE algorithms
# run_type = "bbob"
# ["jde", "jade", "shade"].each{|de_algorithm|
#   out_dir = "results/pal_#{de_algorithm}_n#{pop_size}_#{de_strategy}_#{de_cross}" 

#   if !File.directory?(out_dir)
#     Dir.mkdir(out_dir)
#   end

#   [2, 3, 5, 10, 20, 40].each{|problem_size|
#     max_num_evaluations = problem_size * (10**4)  
#     (1..24).each{|coco_func_num|          
#       instance_out_dir = "#{out_dir}/f#{coco_func_num}_d#{problem_size}"
#       puts instance_out_dir
      
#       if !File.directory?(instance_out_dir)
#         Dir.mkdir(instance_out_dir)
#       end

#       system("python src/createfolders.py #{out_dir}")
#       system("./de -coco_func_num #{coco_func_num} -prosize #{problem_size} -nfe #{max_num_evaluations} -alg #{de_algorithm} -pop #{pop_size} -strategy #{de_strategy} -cross #{de_cross} -sf #{scaling_factor} -cr #{cross_rate} -pop_rate #{pop_rate} -arc_rate #{arc_rate}  -memory_size #{memory_size} -jade_learning_rate #{jade_learning_rate} -p_rate #{p_best_rate}  -tau_sf #{tau_sf} -tau_cr #{tau_cr} -num_cheat_param_sampling #{num_cheat_param_sampling} -out_dir #{instance_out_dir} -run_type #{run_type}")
#     }
#   }
# }
