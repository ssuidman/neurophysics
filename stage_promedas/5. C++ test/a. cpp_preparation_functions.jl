include("../packages.jl")
include("../quickscore_algorithm.jl")

# Data preparation
function prepare_patient_data(case)
    patient_cases_raw = Dict{String,Matrix{Any}}()
    case_nr = parse(Int,case[end])
    file = XLSX.openxlsx("/Users/sam/Large_files/Data/stage_promedas/patient_cases.xlsx")
        for case_nr=1:3
            patient_cases_raw["case $case_nr"] = file["Sheet$case_nr"][1:end, 1:end]
        end
    close(file)
    data_alisa = DataFrame(Float64.(patient_cases_raw[case][8:end,:]),patient_cases_raw[case][6,:])
    idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
    all_sens = Matrix(data_alisa[:,idx_sens])
    pos_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'true'"
    neg_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'false'"
    sens = all_sens[:,pos_idx]
    sensneg = all_sens[:,neg_idx]
    previn = data_alisa[:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]
    return patient_cases_raw, data_alisa, previn, pfmin, pfminneg
end 


function write_cpp_var(previn,pfmin,pfminneg,posterior_BF,dt_exp_sum_log;case_nr=1)
    """
    Write important variables to csv-files for m=7,8,9
    """
    writedlm("variables/cpp/patient404_case_$(case_nr)_previn.csv", previn, ',')
    writedlm("variables/cpp/patient404_case_$(case_nr)_pfmin.csv", pfmin, ',')
    writedlm("variables/cpp/patient404_case_$(case_nr)_pfminneg.csv", pfminneg, ',')
    writedlm("variables/cpp/patient404_case_$(case_nr)_posterior_BF.csv", posterior_BF, ',')
    writedlm("variables/cpp/patient404_case_$(case_nr)_dt_exp_sum_log.csv", dt_exp_sum_log, ',')
end 




