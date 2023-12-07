include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")


# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(18,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data(case)
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick BF thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick MF3 thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick MF4 thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 dict"); 
maximum(abs,posterior .- posterior_)
