include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 3"); (m,n)=size(pfmin); 
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/previn.csv",previn)
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/pfmin.csv",pfmin)
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/pfminneg.csv",pfminneg)
fortran_posterior = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/output_posterior.csv",BigFloat)
fortran_pfplus = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/output_pfplus.csv",BigFloat)
fortran_dt = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/fortran_data/output_time.csv",BigFloat)

method_names = ["prod BF", "prod", "prod Fl128", "prod D64", "prod MF2", "prod MF3", "prod MF4", "prod MF5", "prod BF Fl128","Fortran"]
figure_names_dict = Dict("prod BF"=>"BigFloat","prod"=>"Float64","prod Fl128"=>"Float128","prod D64"=>"Double64","prod MF2"=>"Float64x2","prod MF3"=>"Float64x3","prod MF4"=>"Float64x4","prod MF5"=>"Float64x5","prod BF Fl128"=>"BigFloat128","Fortran"=>"Fortran") 
float_types = [BigFloat, Float64, Float128, Double64, Float64x2, Float64x3, Float64x4, Float64x5, BigFloat, BigFloat]

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
pfplus_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
pfplus_error_dict = Dict(method_name => zeros(m) for method_name in method_names)

for i=1:m
    for method in method_names
        if method!="Fortran"
            pfplus,posterior,dt = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[2,4,5]]; 
            pfplus_dict[method][i],posterior_dict[method][i],dt_dict[method][i] = pfplus[:,1], posterior, dt; 
        else
            posterior_dict[method][i] = fortran_posterior[i,:]
            pfplus_dict[method][i] = fortran_pfplus[i,:]
            dt_dict[method][i] = fortran_dt[i]
        end
        error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["prod BF"][i]) ./ posterior_dict["prod BF"][i])
    end
end
pfplus_BF = zeros(BigFloat,m)
for i=1:m
    pfplus_BF[i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,"trick BF")[2][1,1]
end

# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_time.csv", DataFrame(dt_dict))
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_error.csv", DataFrame(error_dict))
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_posterior.csv", DataFrame(posterior_dict))
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_pfplus_BF.csv", pfplus_BF)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_time.csv",DataFrame)
# posterior_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_posterior.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/high_precision_tools_pfplus_BF.csv")

gr() 
default(fontfamily="Times New Roman") 

p1 = plot()
ylabel!("t (s)") 
for method in method_names[[9,10,3,4,5]]
    plot!(1:m,dt_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10) 
    # plot!(1:m,dt_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10,title="running time") 
end
p2 = plot()
ylabel!(L"$\delta$")
plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
for method in method_names[[10,9,4,5,3]]
    plot!(1:m,error_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10) 
    # plot!(1:m,error_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10,title="relative error") 
end
pa = plot(p1,p2,layout=(1,2),dpi=600) 
plot!(legend = :bottomright,legendfontsize = 7)
xticks!([2,4,6,8], [L"2^2",L"2^{4}",L"2^{6}",L"2^{8}"])
xlabel!(L"#iterations ($2^k$)")
savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/e_high_precision_tools_fig1.png")

p3 = plot()
ylabel!("t (s)") 
for i in [1,8,7,6,5,2]
    method = method_names[i]
    plot!(1:m,dt_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10,color=i) 
    # plot!(1:m,dt_df[:,method],label=figure_names_dict[method],legend=:topleft,yscale=:log10,title="running time",color=i) 
end
p4 = plot()
ylabel!(L"$\delta$")
plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
for i in [2,5,6,7,8]
    method = method_names[i]
    plot!(1:m,error_df[:,method],label=figure_names_dict[method],legend=:topleft,ylim=(10^-95,10^5),yscale=:log10,color=i) 
    # plot!(1:m,error_df[:,method],label=figure_names_dict[method],legend=:topleft,ylim=(10^-95,10^5),yscale=:log10,title="relative error",color=i) 
end
pb = plot(p3,p4,layout=(1,2),dpi=600) 
plot!(legend = :bottomright,legendfontsize = 7)
xticks!([2,4,6,8], [L"2^2",L"2^{4}",L"2^{6}",L"2^{8}"])
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/e_high_precision_tools_fig2.png")


p5 = plot()
ylabel!("t (s)") 
for i in [1,8,7,3,6,4,5,2]
    method = method_names[i]
    plot!(1:m,dt_df[:,method],label=figure_names_dict[method],line=(1,:solid),yscale=:log10,color=i) 
    # plot!(1:m,dt_df[:,method],label=figure_names_dict[method],line=(1,:solid),yscale=:log10,title="running time",color=i) 
    a,b = linear_fit(5:9,log10.(dt_df[5:9,method]))
    plot!(5:m+22,10 .^(a .+ b.*[5:m+22...]),label=nothing,line=(1,:dash),yscale=:log10,color=i) 
end
p6 = plot()
ylabel!(L"$\delta$")
plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",line=(1,:dashdot),legend=:topleft,yscale=:log10,linestyle=:dash,color="black") 
# plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",line=(1,:dashdot),legend=:topleft,yscale=:log10,title="relative error",linestyle=:dash,color="black") 
for i in [2,4,5,3,6,7,8]
    method = method_names[i]
    plot!(1:m,error_df[:,method],label=figure_names_dict[method],line=(1,:solid),legend=:topleft,ylim=(10^-95,10^10),yscale=:log10,color=i) 
    # plot!(1:m,error_df[:,method],label=figure_names_dict[method],line=(1,:solid),legend=:topleft,ylim=(10^-95,10^10),yscale=:log10,title="relative error",color=i) 
    a,b = linear_fit(5:9,log10.(error_df[5:9,method]))
    plot!(5:m+22,10 .^(a .+ b.*[5:m+22...]),label=nothing,line=(1,:dash),yscale=:log10,color=i)   
end
a,b = linear_fit(5:9,log10.(pfplus_BF[5:9]))
plot!(5:m+22,10 .^(a .+ b.*[5:m+22...]),label=nothing,line=(1,:dash),yscale=:log10,color="grey")
pc = plot(p5,p6,layout=(1,2),dpi=600) 
plot!(legend = :bottomright,legendfontsize = 7)
xticks!([5,10,15,20,25,30], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}",L"2^{25}",L"2^{30}"])
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/e_high_precision_tools_fig3.png")







