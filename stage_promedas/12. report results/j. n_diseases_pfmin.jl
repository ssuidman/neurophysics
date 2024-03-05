include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

# Hier begin analyse waarin ik laat zien dat hoe kleiner het aantal disease n, hoe kleiner pfmin 
previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick MF5 thread");  

pfmin_test = copy(pfmin)
pfmin_test[pfmin.==1] = 1 .- rand(sum(pfmin.==1)).*0.0001
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/n_diseases_pfmin.csv",pfmin_test)
pfmin_test = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/n_diseases_pfmin.csv")


pfplus_dict = Dict{String,MultiFloat{Float64, 4}}()
i_values = [3,5,8,10,20,40,60,80,100:25:1000...,n]
for i in i_values
    println("i=$i")
    pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn[1:i], pfmin_test[:,1:i], pfminneg[1:i],"trick MF4 thread"); 
    pfplus_dict["$i"] = pfplus[1]
end
# pfplus_df = DataFrame(pfplus_dict)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/n_diseases_pfplus.csv",pfplus_df)
pfplus_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/n_diseases_pfplus.csv",DataFrame)

pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick MF4 thread"); 

p = scatter([n],[pfplus[1]],yscale=:log10,ylims=(10^-60,10^0),yticks=[10^-60,10^-40,10^-20,10^0],label="Promedas",color="blue",markersize=5,legend=:topleft,dpi=600)
plot!(i_values,Vector(pfplus_df[1,string.([i_values...])]),label="random",yscale=:log10,color="red",xlabel="First i diseases",ylabel=latexstring("p(F^+,F^-)"))
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/j_n_diseases_pfmin.png")
