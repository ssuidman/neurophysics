include("../packages.jl")
include("../quickscore_algorithm.jl")
include("a. cpp_preparation_functions.jl")

for case in ["case 1","case 2","case 3"]
    patient_cases_raw, data_alisa, previn, pfmin, pfminneg = prepare_patient_data(case)
    pfplus_BF, P_joint_BF, posterior_BF, dt_BF = quickscore(previn, pfmin, pfminneg, "prod BF") 
    pfplus_exp_sum_log, P_joint_exp_sum_log, posterior_exp_sum_log, dt_exp_sum_log = quickscore(previn, pfmin, pfminneg, "exp-sum-log") 
    # DEZE REGEL VERVANGT DE EXCEL BESTANDEN DUS OPPASSEN!!!
    # write_cpp_var(previn,pfmin,pfminneg,posterior_BF,dt_exp_sum_log,case_nr=parse(Int,case[end]))
end