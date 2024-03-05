include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 3"); (m,n)=size(pfmin); 

method_names = ["prod BF", "trick BF", "prod MF2", "trick MF2", "prod", "trick"]
float_types = [BigFloat, BigFloat, Float64x2, Float64x2, Float64, Float64]
names_dict = Dict("prod BF"=>"BigFloat","trick BF"=>"BigFloat","prod MF2"=>"Float64x2","trick MF2"=>"Float64x2","prod"=>"Float64","trick"=>"Float64") 
color_dict = Dict("prod BF"=>1,"trick BF"=>1,"prod MF2"=>2,"trick MF2"=>2,"prod"=>3,"trick"=>3) 

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
for i=1:m
    for method in method_names
        posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
        error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["prod BF"][i]) ./ posterior_dict["prod BF"][i])
    end
end
pfplus_BF = zeros(BigFloat,m)
for i=1:m
    pfplus_BF[i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,"trick BF")[2][1,1]
end

# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_time.csv", DataFrame(dt_dict))
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_error.csv", DataFrame(error_dict))
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_pfplus_BF.csv",pfplus_BF)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_time.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/trick_pfplus_BF.csv")

gr() 
default(fontfamily="Times New Roman") 

p1 = plot()
ylabel!("t (s)") 
for method in method_names
    if occursin("trick",method) 
        plot!(1:m,dt_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-3.95,10^3),
            yscale=:log10,
            # title="running time",
            color=color_dict[method]
            ) 
    else
        plot!(1:m,dt_df[:,method],
            label=false,
            legend=:topleft,
            yscale=:log10,
            # title="running time",
            color=color_dict[method],
            linestyle=:dash
        ) 
    end
end
p2 = plot()
ylabel!(L"$\delta$")
plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:m,pfplus_BF,label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
for method in method_names[[2,3,4,5,6]]
    if occursin("trick",method)
        plot!(1:m,error_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
            ) 
        else
            plot!(1:m,error_df[:,method],
            label=false,
            legend=:topleft,
            ylim=(10^-90,10^5),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method],
            linestyle=:dash
        ) 
    end
end
pa = plot(p1,p2,layout=(1,2),dpi=600) 
plot!(legend = :bottomright,legendfontsize = 7)
xticks!([2,4,6,8], [L"2^2",L"2^{4}",L"2^{6}",L"2^{8}"])
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/f_trick.png")

