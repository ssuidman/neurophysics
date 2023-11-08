include("../packages.jl")
include("../quickscore_algorithm.jl")
include("../run_one_time.jl")

# This file is used such that "c. helpful_tools.jl" can use it to look at precompiling 
previn, pfmin, pfminneg, pfplus_, P_joint_, posterior_, dt_, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3);
pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"exp-sum-log");