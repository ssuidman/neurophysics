include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../0. useful/run_one_time.jl")

# This file is used such that "c. helpful_tools.jl" can use it to look at precompiling 
previn, pfmin, pfminneg, pfplus_, P_joint_, posterior_, dt_, prev, prevminneg, myset = run_one_time_var(m=20,n_myset=3);
pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"exp-sum-log");
pfplus_matrix_128, pfplus_128, P_joint_128, posterior_128, dt_128 = quickscore(previn,pfmin,pfminneg,"Fl128 trick thread");
pfplus_matrix_MF2, pfplus_MF2, P_joint_MF2, posterior_MF2, dt_MF2 = quickscore(previn,pfmin,pfminneg,"MF2 trick thread");
pfplus_matrix_MF3, pfplus_MF3, P_joint_MF3, posterior_MF3, dt_MF3 = quickscore(previn,pfmin,pfminneg,"MF3 trick thread");
pfplus_matrix_MF4, pfplus_MF4, P_joint_MF4, posterior_MF4, dt_MF4 = quickscore(previn,pfmin,pfminneg,"MF4 trick thread");
pfplus_matrix_BF, pfplus_BF, P_joint_BF, posterior_BF, dt_BF = quickscore(previn,pfmin,pfminneg,"BF trick thread");


maximum(abs,pfplus_BF.-pfplus_128)
maximum(abs,pfplus_BF.-pfplus_MF2)
maximum(abs,pfplus_BF.-pfplus_MF3)
maximum(abs,pfplus_BF.-pfplus_MF4)

maximum(abs,(posterior_BF.-posterior_128)./posterior_BF)
maximum(abs,(posterior_BF.-posterior_MF2)./posterior_BF)
maximum(abs,(posterior_BF.-posterior_MF3)./posterior_BF)
maximum(abs,(posterior_BF.-posterior_MF4)./posterior_BF)