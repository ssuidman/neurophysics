include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")
include("0. quickscore_algorithm_noisy_max_report.jl")

patientprev,sens_normal,sens_medium,sens_extreme = prepare_patient_data("case 8"); n = size(patientprev,1); m1=size(sens_extreme,2); m2=size(sens_medium,2); m=m1+m2; 

method_names = ["trick BF thread","trick thread","trick MF2 thread","trick MF3 thread","trick MF4 thread","trick MF5 thread"]
float_types = [BigFloat,Float64,Float64x2,Float64x3,Float64x4,Float64x5]
names_dict = Dict("trick BF thread"=>"BigFloat","trick thread"=>"Float64","trick MF2 thread"=>"Float64x2","trick MF3 thread"=>"Float64x3","trick MF4 thread"=>"Float64x4","trick MF5 thread"=>"Float64x5") 
color_dict = Dict("trick BF thread"=>1,"trick thread"=>2,"trick MF2 thread"=>3,"trick MF3 thread"=>4,"trick MF4 thread"=>5,"trick MF5 thread"=>6) 

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
pfplus_BF = zeros(BigFloat,m)

for i=1:m1
    for method in method_names
        # quickscore_noisy_max(patientprev, sens_normal, sens_medium[:,[],:], sens_extreme[:,1:i], method)[[3,4]]
        pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_noisy_max(patientprev, sens_normal, sens_medium[:,[],:], sens_extreme[:,1:i], method)[[1,3,4]]
        error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
        pfplus_BF[i] = pfplus[1,1]
    end
end
for i=1:m2
    for method in method_names
        pfplus,posterior_dict[method][i+m1],dt_dict[method][i+m1] = quickscore_noisy_max(patientprev, sens_normal, sens_medium[:,1:i,:], sens_extreme, method)[[2,3,4]]
        error_dict[method][i+m1] = maximum(abs,(posterior_dict[method][i+m1] .- posterior_dict["trick BF thread"][i+m1]) ./ posterior_dict["trick BF thread"][i+m1])
        pfplus_BF[i+m1] = pfplus[1,1]
    end
end
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_time.csv",DataFrame(dt_dict))
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_error.csv", DataFrame(error_dict))
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_pfplus_BF.csv",pfplus_BF)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_time.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_noisy_max_pfplus_BF.csv")

gr() 
default(fontfamily="Times New Roman") 

p1 = plot()
ylabel!("t (s)") 
for method in method_names
    println(method)
    if occursin("thread",method) 
        plot!(1:m,dt_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-3.75,10^3),
            yscale=:log10,
            # title="running time",
            color=color_dict[method]
        ) 
    else
        # plot!(1:m,dt_df[:,method],
        #     label=false,
        #     legend=:topleft,
        #     yscale=:log10,
        #     title="running time",
        #     color=color_dict[method],
        #     linestyle=:dash
        # ) 
    end
end
xticks!([5,10,15,18], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{18}"])
plot!(legend = :bottomright,legendfontsize = 7)
p2 = plot()
ylabel!(L"$\delta$")
plot!(1:18,pfplus_BF[1:18],label=L"p(F^{++},F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:18,pfplus_BF[1:18],label=L"p(F^{++},F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
for i in 2:6
    method = method_names[i]
    if occursin("thread",method) & !occursin("BF",method)
        hline!([eps(float_types[i])],color=color_dict[method],linestyle=:dashdot,label=false)
        annotate!(15, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(names_dict[method])}}}"), :left, :bottom, 7, :black))
        plot!(1:18,error_df[1:18,method],
            label=names_dict[method],
            legend=:topleft,
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
        ) 
    elseif method=="trick BF thread"
        plot!(3:18,error_df[3:18,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-100,10^2),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
        ) 
    else
        plot!(1:18,error_df[1:18,method],
            label=false,
            legend=:topleft,
            ylim=(10^-100,10^2),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method],
            line=(2,:dash)
        ) 
    end
end
xticks!([5,10,15,18], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{18}"])
plot!(legend = :bottom,legendfontsize = 6)
pa = plot(p1,p2,layout=(1,2),dpi=600) 
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/h_thread_noisy_max.png")
