include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin); 

method_names = ["trick BF", "trick dict BF", "trick matrix BF", "trick MF2", "trick dict MF2", "trick", "trick dict"]
float_types = [BigFloat, BigFloat, BigFloat, Float64x2, Float64x2, Float64, Float64]
# names_dict = Dict("prod BF"=>"BigFloat","prod dict BF"=>"BigFloat","prod MF2"=>"Float64x2","prod dict MF2"=>"Float64x2","prod"=>"Float64","prod dict"=>"Float64") 
names_dict = Dict("trick BF"=>"BigFloat","trick dict BF"=>"BigFloat","trick matrix BF"=>"BigFloat","trick MF2"=>"Float64x2","trick dict MF2"=>"Float64x2","trick"=>"Float64","trick dict"=>"Float64") 
color_dict = Dict("trick BF"=>1,"trick dict BF"=>1,"trick matrix BF"=>1,"trick MF2"=>2,"trick dict MF2"=>2,"trick"=>3,"trick dict"=>3) 

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
# First without the "dict" method:
for i=1:m
    for method in method_names[[1,3,4,5,6,7]]
        println(method)
        posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
        error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF"][i]) ./ posterior_dict["trick BF"][i])
    end
end
# Then add the dict method up to m-1
for i=1:m-1
    posterior_dict["trick dict BF"][i],dt_dict["trick dict BF"][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,"trick dict BF")[[4,5]]; 
    error_dict["trick dict BF"][i] = maximum(abs,(posterior_dict["trick dict BF"][i] .- posterior_dict["trick BF"][i]) ./ posterior_dict["trick BF"][i])
end
# Then add the last term that blows up manually 
posterior_dict["trick dict BF"][m],dt_dict["trick dict BF"][m] = zeros(n), 9233.0;
error_dict["trick dict BF"][m] = error_dict["trick matrix BF"][m]

pfplus_BF = zeros(BigFloat,m)
for i=1:m
    pfplus_BF[i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,"trick BF thread")[2][1,1]
end
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_time.csv", DataFrame(dt_dict))
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_error.csv", DataFrame(error_dict))
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_pfplus_BF.csv",pfplus_BF)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_time.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/dict_pfplus_BF.csv")

gr() 
default(fontfamily="Times New Roman") 

p1 = plot()
ylabel!("t (s)") 
for method in method_names[[1,2,4,5,6,7]]
    println(method)
    if occursin("dict",method) 
        plot!(1:m,dt_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-4.75,10^3),
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
for method in method_names[[2,4,5,6,7]]
    if occursin("dict",method) & !occursin("BF",method)
        plot!(1:m,error_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
        ) 
    elseif occursin("BF",method)
        plot!(11:m,error_df[11:m,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-70,10^2),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method],
        ) 
    else
        plot!(1:m,error_df[:,method],
            label=false,
            legend=:topleft,
            ylim=(10^-70,10^2),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method],
            line=(1,:dash)
        ) 
    end
end
pa = plot(p1,p2,layout=(1,2),dpi=600) 
plot!(legend = :bottomleft,legendfontsize = 7)
xticks!([2,4,6,8,10,12,14], [L"2^{2}",L"2^{4}",L"2^{6}",L"2^{8}",L"2^{10}",L"2^{12}",L"2^{14}"])
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/g_dict.png")

