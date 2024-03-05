include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 4"); (m,n)=size(pfmin); 

pfplus_matrix_32, pfplus_32, P_joint_32, posterior_32, dt_32, dt_array_32 = quickscore(previn, pfmin, pfminneg,"trick Fl32 linear");  
pfplus_matrix_64, pfplus_64, P_joint_64, posterior_64, dt_64, dt_array_64 = quickscore(previn, pfmin, pfminneg,"trick linear");  
pfplus_matrix_128, pfplus_128, P_joint_128, posterior_128, dt_128, dt_array_128 = quickscore(previn, pfmin, pfminneg,"trick Fl128 linear");  
pfplus_matrix_BF, pfplus_BF, P_joint_BF, posterior_BF, dt_BF, dt_array_BF = quickscore(previn, pfmin, pfminneg,"trick BF linear");  

cumsum_32 = cumsum(pfplus_matrix_32[:,1,1])[2 .^[0:13...]]
cumsum_64 = cumsum(pfplus_matrix_64[:,1,1])[2 .^[0:13...]]
cumsum_128 = cumsum(pfplus_matrix_128[:,1,1])[2 .^[0:13...]]
cumsum_BF = cumsum(pfplus_matrix_BF[:,1,1])[2 .^[0:13...]]

# First 2^0, 2^1, 2^2, 2^3, 2^4 terms of numerical evalution of quickscore
sum(pfplus_matrix_BF[1:2^0,1,1])
sum(pfplus_matrix_BF[1:2^1,1,1])
sum(pfplus_matrix_BF[1:2^2,1,1])
sum(pfplus_matrix_BF[1:2^3,1,1])
sum(pfplus_matrix_BF[1:2^4,1,1])

# Table 2^k
relative_error_table = Float64.(hcat(0:13,cumsum_64,cumsum_BF,abs.(cumsum_32.-cumsum_BF)./cumsum_BF,abs.(cumsum_64.-cumsum_BF)./cumsum_BF,abs.(cumsum_128.-cumsum_BF)./cumsum_BF))
# writedlm("variables/report_results/relative_error_table.csv",relative_error_table,",")

