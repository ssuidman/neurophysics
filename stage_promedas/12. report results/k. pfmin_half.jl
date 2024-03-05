include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

# Hier begin analyse waarin ik laat zien dat hoe kleiner het aantal disease n, hoe kleiner pfmin 
previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin[1:15,:], pfminneg,"trick MF3 thread");  

# rand_data = (rand(size(pfmin)...).*2 .-1).*0.05
# rand_data[pfmin.==1] .= 0 
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_random_data.csv",rand_data)
rand_data = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_random_data.csv")

k_values = 0.05:0.05:0.95
# k_values = [0.1,0.3,0.5,0.7,0.9]
pfplus_test_dict = Dict{String,MultiFloat{Float64, 4}}()
pfplus_rand_dict = Dict{String,MultiFloat{Float64, 4}}()
for k in k_values
    println("k = $k")
    pfmin_test, pfmin_rand = copy(pfmin), copy(pfmin)
    pfmin_test .= k
    pfmin_rand .= k .+ rand_data
    pfmin_test[pfmin.==1] .= 1
    pfmin_rand[pfmin.==1] .= 1
    # pretty_table(pfmin_test)
    # pretty_table(pfmin_rand)
    pfplus_matrix, pfplus_test, P_joint, posterior, dt = quickscore_report(previn, pfmin_test, pfminneg,"trick MF3 thread");  
    pfplus_matrix, pfplus_rand, P_joint, posterior, dt = quickscore_report(previn, pfmin_rand, pfminneg,"trick MF3 thread");  
    pfplus_test_dict["$k"] = pfplus_test[1]
    pfplus_rand_dict["$k"] = pfplus_rand[1]
end
# pfplus_test_df = DataFrame(pfplus_test_dict)
# pfplus_rand_df = DataFrame(pfplus_rand_dict)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_pfplus_test.csv",pfplus_test_df)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_pfplus_rand.csv",pfplus_rand_df)
pfplus_test_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_pfplus_test.csv",DataFrame)
pfplus_rand_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/pfmin_half_pfplus_rand.csv",DataFrame)

pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick MF3 thread"); 
p = scatter([mean(1 .- pfmin[pfmin.!=1])],[pfplus[1]],yscale=:log10,label="Promedas mean",xlabel=L"\tilde{w}_{ji}",ylabel=L"p(F^+,F^-)",ylims=(10^-31,1),yticks=[10^-30,10^-20,10^-10,10^0],color="blue",markersize=5,legend=:topleft,dpi=600)
plot!(1 .- k_values,Vector(pfplus_test_df[1,string.([k_values...])]),color="red",linestyle=:solid,label="constant")
plot!(1 .- k_values,Vector(pfplus_rand_df[1,string.([k_values...])]),color="green",linestyle=:dash,label="constant + noise")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/k_pfmin_half.png")



