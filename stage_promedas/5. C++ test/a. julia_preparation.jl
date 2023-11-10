include("../packages.jl")
include("../quickscore_algorithm.jl")

# Data preparation
function prepare_cpp_var()
    matrix_alisa = Dict{String,Matrix{Any}}()
    file = XLSX.openxlsx("/Users/sam/Large_files/Data/stage_promedas/patient_cases.xlsx")
        matrix_alisa["case 1"] = file["Sheet1"][1:end, 1:end]
    close(file)
    case = "case 1"
    data_alisa = DataFrame(Float64.(matrix_alisa[case][8:end,:]),matrix_alisa[case][6,:])
    idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
    all_sens = Matrix(data_alisa[:,idx_sens])
    pos_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'true'"
    neg_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'false'"
    sens = all_sens[:,pos_idx]
    sensneg = all_sens[:,neg_idx]
    previn = data_alisa[:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]
    return previn, pfmin, pfminneg 
end 

previn, pfmin, pfminneg = prepare_cpp_var()
pfplus_BF, P_joint_BF, posterior_BF, dt_BF = quickscore(previn, pfmin, pfminneg, "prod BF") 
pfplus_exp_sum_log, P_joint_exp_sum_log, posterior_exp_sum_log, dt_exp_sum_log = quickscore(previn, pfmin, pfminneg, "exp-sum-log") 

function write_cpp_var(previn,pfmin,pfminneg,posterior_BF,dt_exp_sum_log)
    # Write important variables to csv-files 
    writedlm("variables/cpp_preparation/patient404_case_1_previn.csv", previn, ',')
    writedlm("variables/cpp_preparation/patient404_case_1_pfmin.csv", pfmin, ',')
    writedlm("variables/cpp_preparation/patient404_case_1_pfminneg.csv", pfminneg, ',')
    writedlm("variables/cpp_preparation/patient404_case_1_posterior_BF.csv", posterior_BF, ',')
    writedlm("variables/cpp_preparation/patient404_case_1_dt_exp_sum_log.csv", dt_exp_sum_log, ',')
end 
write_cpp_var(previn,pfmin,pfminneg,posterior_BF,dt_exp_sum_log)



