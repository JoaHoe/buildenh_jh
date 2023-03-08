##name of script: func_loadLib_jh.R
cat("version_number= ",v_nr,"\n")
#purpose: collection of all functions for package 'buildenh'
#author: Joachim Höhle
#GNU General Public License (GPL)

loadLib_jh <- function() { 
  source("./func/func_adjust_coord.R")
  source("./func/func_corner_name.R")
  source("./func/func_design_mat.R")
  source("./func/func_design_mat2.R")
  source("./func/func_det_of_angle.R")
  source("./func/func_detect_meas1.R")
  source("./func/func_dist_PC.R")
  source("./func/func_dist_v.R")
  source("./func/func_dist_v2.R")
  source("./func/func_fe.R")
  source("./func/func_line_reduce.R")
  source("./func/func_locator2.R")
  source("./func/func_mean_line.R")
  source("./func/func_new_centre.R")
  source("./func/func_new_centre_auto.R")
  source("./func/func_PC_segment_4.R")
  source("./func/func_rectang_lines.R")
  source("./func/func_reduce_pointset.R")
  source("./func/func_test_gaps_onedirection.R")
  source("./func/func_detect_thet_len.R")
  source("./func/func_trans_ortho.R")
  source("./func/func_trans_H_res.R")
  source("./func/func_w_av.R")
} #end of function 'loadLib_jh()'

#end of script 'func_loadLib_jh.R'