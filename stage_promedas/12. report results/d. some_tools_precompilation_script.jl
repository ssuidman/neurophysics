include("0. quickscore_algorithm_report.jl")

function prepare_patient_data(case)
    patient_cases_raw = Dict{String,Matrix{Any}}()
    case_nr = parse(Int,case[end])
    file = XLSX.openxlsx("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_cases.xlsx")
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

    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    println("n = $n \t m+ = $m \t m- = $(size(sensneg,2))")
    return previn, pfmin, pfminneg
end 


# This file is used such that "c. helpful_tools.jl" can use it to look at precompiling 
previn, pfmin, pfminneg = prepare_patient_data("case 3")
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"prod");  