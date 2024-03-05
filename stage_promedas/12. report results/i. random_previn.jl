include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")


######################################################################################################
##################################### figuur maken  ##################################################
######################################################################################################
previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick MF5 thread");  
k_values = -1:1:2
previn_dict = Dict{String,Vector{Float64}}()
for k in k_values
    println(k)
    log_x = log10.(previn) 
    u,o = (mean(log_x).+k)./-maximum(abs.(log_x)), std(log_x./-maximum(abs.(log_x)))
    alpha = ((1-u)/o^2-1/u)*u^2
    beta = alpha*(1/u-1)
    previn_test = 10 .^(rand(Beta(alpha,beta),size(previn)...).*-maximum(abs.(log_x)))
    previn_dict["$k"] = previn_test
end 
# previn_df = DataFrame(previn_dict)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_previn.csv",previn_df)
previn_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_previn.csv",DataFrame)

h_vec = Vector{Any}(undef,4)
for (i,k) in enumerate(k_values)
    previn_test = previn_df[:,"$k"]
    h = histogram(log10.(previn),alpha=0.7,normalize=true,bins=30,label=false,color="blue",xlim=(-8,0))
    if sign(k)==-1; nr="-$(abs(k))"; elseif(sign(k))==1; nr="+$(abs(k))"; else; nr="+$(abs(k))"; end; 
    histogram!(log10.(previn_test),alpha=0.4,normalize=true,bins=30,label=false,color="red")
    annotate!(-4,0.8,text(latexstring("\\mu=$(Int(round(mean(log10.(previn)),sigdigits=1)))"),:center,10,color=:blue))
    if k<=0
        annotate!(-6,0.45,text(latexstring("\\mu$(nr)"),:left,10,color=:red))
    else
        annotate!(-3,0.45,text(latexstring("\\mu$(nr)"),:left,10,color=:red))
    end
    if i>2; histogram!(xlabel=latexstring("\\log_{10}\\left[p(d_i)\\right]")); end; 
    histogram!(xlabel=latexstring("\\log_{10}\\left[p(d_i)\\right]")); 
    h_vec[i] = h
end
p1 = plot(h_vec...,layout=(2,2),tick_direction=:in,ylims=(0,1),dpi=600)
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/i_random_previn_fig1.png")

######################################################################################################
##################################### data plotten  ##################################################
######################################################################################################

k_values2 = -3:0.25:4
previn_dict2 = Dict{String,Vector{Float64}}()
for k in k_values2
    println("k=$k")
    log_x = log10.(previn) 
    u,o = (mean(log_x).+k)./-maximum(abs.(log_x)), std(log_x./-maximum(abs.(log_x)))
    alpha = ((1-u)/o^2-1/u)*u^2
    beta = alpha*(1/u-1)
    previn_test = 10 .^(rand(Beta(alpha,beta),size(previn)...).*-maximum(abs.(log_x)))
    previn_dict2["$k"] = previn_test
end
# previn_df2 = DataFrame(previn_dict2)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_previn2.csv",previn_df2)
previn_df2 = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_previn2.csv",DataFrame)

pfplus_dict = Dict{String,MultiFloat{Float64, 3}}()
for k in k_values2
    println("k=$k")
    previn_test = previn_df2[:,"$k"]
    pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn_test, pfmin, pfminneg,"trick MF3 thread");  
    pfplus_dict["$k"] = pfplus[1]
end 
# pfplus_df = DataFrame(pfplus_dict)
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_pfplus.csv",pfplus_df)
pfplus_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/random_pfplus.csv",DataFrame)

pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick MF3 thread");  

p2 = scatter([mean(log10.(previn))],[pfplus[1]],yscale=:log10,ylims=(10^-35,10^0),yticks=[10^-30,10^-20,10^-10,10^0],label="Promedas",color="blue",markersize=5,xlims=(-8,0),legend=:topleft,dpi=600)
plot!(mean(log10.(previn)).+k_values2,Vector(pfplus_df[1,string.([k_values2...])]),label="random",yscale=:log10,color="red",xlabel=L"\overline{\log_{10}\left[p(d_i)\right]}",ylabel=latexstring("p(F^+,F^-)"))
# plot!(mean(log10.(previn)).+k_values2,Vector(pfplus_df[1,string.([k_values2...])]),label="random",yscale=:log10,color="red",xlabel=latexstring("\\left\\langle \\log_{10}\\left[p(d_i)\\right] \\right\\rangle"),ylabel=latexstring("p(F^+,F^-)"))
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/i_random_previn_fig2.png")


