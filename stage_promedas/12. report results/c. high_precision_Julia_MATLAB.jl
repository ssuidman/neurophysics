include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 3"); 
# dt_array_Float64 = zeros(Float64,9)
# dt_array_BF = zeros(BigFloat,9)
# for i=1:9
#     pfplus_matrix, pfplus, P_joint, posterior, dt, dt_array_not = quickscore_report(previn,pfmin[1:i,:],pfminneg,"prod")
#     dt_array_Float64[i] = dt
#     pfplus_matrix, pfplus, P_joint, posterior, dt, dt_array_not = quickscore_report(previn,pfmin[1:i,:],pfminneg,"prod BF")
#     dt_array_BF[i] = dt
# end

# To not have to run the data again the running times are shown here 
dt_array_Float64 = [0.038284917, 0.007979291, 0.027357709, 0.054937916, 0.050511458, 0.142536708, 0.207391584, 0.352631959, 0.618278208] 
dt_array_BF = [1.030957, 0.95503025, 1.762607333, 3.919750125, 8.216613417, 19.829208291, 44.959633959, 61.638607083, 118.417239833] 
dt_array_matlab_double = [0.0251, 0.0430, 0.0802, 0.1386, 0.2443, 0.4569, 0.8685, 1.7031, 3.3893] 
dt_array_matlab_Mutliprecision = [4.20301570, 7.95905320, 15.5147344, 30.8469886, 69.9844129, 139.7692136, 263.7503738, 522.1046762, 1021.7811609]
# # Save the running times in a csv-file in variables-folder 
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_julia_matlab.csv",hcat(dt_array_Float64,dt_array_BF,dt_array_matlab_double,dt_array_matlab_Mutliprecision),',')

dt_matrix = DataFrame(readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_julia_matlab.csv",','),["Float64","BigFloat","MATLAB double","MATLAB Multiprecision"])

gr() 
default(fontfamily="Times New Roman") 
p1 = plot(1:9,dt_matrix[:,"Float64"],ylims = (10^-2.5,10^3.2),label="double precision",legend=:topleft,yscale=:log10) 
# p1 = plot(1:9,dt_matrix[:,"Float64"],ylims = (10^-2.5,10^3.2),label="double precision",legend=:topleft,yscale=:log10,title="Julia") 
plot!(1:9,dt_matrix[:,"BigFloat"],label="arbitrary precision",yscale=:log10) 
p2 = plot(1:9,dt_matrix[:,"MATLAB double"],label="double precision",legend=:topleft,yscale=:log10) 
# p2 = plot(1:9,dt_matrix[:,"MATLAB double"],label="double precision",legend=:topleft,yscale=:log10,title="MATLAB") 
plot!(1:9,dt_matrix[:,"MATLAB Multiprecision"],ylims = (10^-2.5,10^3.2),label="arbitrary precision",yscale=:log10) 
p = plot(p1,p2,layout=(1,2),dpi=600) 
xticks!([2,4,6,8], [L"2^2",L"2^{4}",L"2^{6}",L"2^{8}"])
xlabel!(L"#iterations ($2^k$)")
# xlabel!(L"positive findings ($m^+$)") 
ylabel!("t (s)") 
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/c_high_precision_plot.png")