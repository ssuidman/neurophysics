include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn, pfmin, pfminneg = prepare_patient_data("case 6"); (m,n)=size(pfmin);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin[1:18,:], pfminneg,"trick MF5 thread");  


method_names = ["trick BF thread", "trick MF2", "trick MF2 thread", "trick MF3", "trick MF3 thread", "trick MF4 thread","trick MF5 thread", "trick", "trick thread"]
float_types = [BigFloat, Float64x2, Float64x2, Float64x3, Float64x3, Float64x4, Float64x5, Float64, Float64]
names_dict = Dict("trick BF thread"=>"BigFloat","trick MF2"=>"Float64x2","trick MF2 thread"=>"Float64x2","trick MF3"=>"Float64x3","trick MF3 thread"=>"Float64x3","trick MF4 thread"=>"Float64x4","trick MF5 thread"=>"Float64x5","trick"=>"Float64","trick thread"=>"Float64") 
eps_names_dict = Dict("trick BF thread"=>"BF","trick MF2"=>"64x2","trick MF2 thread"=>"64x2","trick MF3"=>"64x3","trick MF3 thread"=>"64x3","trick MF4 thread"=>"64x4","trick MF5 thread"=>"64x5","trick"=>"64","trick thread"=>"64") 
color_dict = Dict("trick BF thread"=>1,"trick MF2"=>2,"trick MF2 thread"=>2,"trick MF3"=>3,"trick MF3 thread"=>3,"trick MF4 thread"=>4,"trick MF5 thread"=>5,"trick"=>6,"trick thread"=>6) 

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
pfplus_BF = zeros(BigFloat,18)
# For the BF methods only up to m=18
for i=1:18
    for method in method_names[[1,7]]
        if occursin("BF",method) && i>18
            posterior_dict[method][i],dt_dict[method][i] = zeros(n), 0.0 
        elseif occursin("BF",method) && i<=18
            pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[2,4,5]]; 
            error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
            pfplus_BF[i] = pfplus[1,1]
        elseif !occursin("BF",method) && i<=18
            posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
            error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
        else
            posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
        end
    end
end
# Then calculate BigFloat pfplus 
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_time.csv",DataFrame(dt_dict)[:,["trick MF5 thread"]])
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_error.csv", DataFrame(error_dict)[:,["trick MF5 thread"]])
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_pfplus_BF.csv",pfplus_BF)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_time.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_pfplus_BF.csv")

gr() 
default(fontfamily="Times New Roman") 

p1 = plot()
ylabel!("t (s)") 
for method in method_names[[1,7,6,5,3,9]]
    if occursin("thread",method) && any(occursin.(["MF5","BF"],method))
        plot!(1:18,dt_df[1:18,method],
        label=names_dict[method],
        legend=:topleft,
        ylim=(10^-3.75,10^3),
        yscale=:log10,
        # title="running time",
        color=color_dict[method]
        ) 
    else 
        println(method," ",dt_df[17:22,method]," 2")
        plot!(1:m,dt_df[:,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-3.75,10^3),
            yscale=:log10,
            # title="running time",
            color=color_dict[method]
        ) 
    end
end
xticks!([5,10,15,20], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}"])
plot!(legend = :bottomright,legendfontsize = 7)
p2 = plot()
ylabel!(L"$\delta$")
plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
for i in [8,9,2,3,4,5,6,7,1]
    method = method_names[i]
    if occursin("thread",method) & !occursin("BF",method)
        hline!([eps(float_types[i])],color=color_dict[method],linestyle=:dashdot,label=false)
        annotate!(16.5, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
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
    end
end
xticks!([5,10,15,18], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{18}"])
plot!(legend = :bottom,legendfontsize = 6)
pa = plot(p1,p2,layout=(1,2),dpi=600) 
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/h_thread_fig1.png")



p3 = plot()
ylabel!("t (s)") 
for method in method_names[[1,7,6,4,5,2,3,8,9]]
    if occursin("thread",method) 
        if occursin("BF",method) || occursin("MF5",method); m=18; else; m=22; end; println(method," ",m)
        a,b = linear_fit(10:m,log10.(dt_df[10:m,method]))
        plot!(10:35,10 .^(a .+ b.*[10:35...]),label=nothing,line=(1,:dash),yscale=:log10,color=color_dict[method])
        println(1:m," ",dt_df[1:m,method])
        plot!(1:m,dt_df[1:m,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-4,10^8),
            yscale=:log10,
            # title="running time",
            color=color_dict[method]
        ) 
    end
end
xticks!([5,10,15,20,25,30,35], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}",L"2^{25}",L"2^{30}",L"2^{35}"])
plot!(legend = :bottomright,legendfontsize = 7)
p4 = plot()
ylabel!(L"$\delta$")
plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
a,b = linear_fit(8:18,log10.(pfplus_BF[8:18]))
plot!(8:35,10 .^(a .+ b.*[8:35...]),label=nothing,line=(1,:dash),yscale=:log10,color="black")
for i in [8,9,2,3,4,5,6,7,1]
    method = method_names[i]
    if occursin("thread",method) & !occursin("BF",method)
        hline!([eps(float_types[i])],color=color_dict[method],linestyle=:dashdot,label=false)
        annotate!(32.5, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
        plot!(1:18,error_df[1:18,method],
            label=names_dict[method],
            legend=:topleft,
            yscale=:log10,
            ylim=(10^-80,10^2),
            # title="relative error",
            color=color_dict[method]
        ) 
    # elseif method=="trick BF thread"
    #     setprecision(BigFloat,300) # because this is also done in quickscore algorithm
    #     hline!([eps(float_types[i])],color=color_dict[method],linestyle=:dashdot,label=false)
    #     annotate!(32.5, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
    #     plot!(3:18,error_df[3:18,method],
    #         label=names_dict[method],
    #         legend=:topleft,
    #         ylim=(10^-110,10^2),
    #         yscale=:log10,
    #         # title="relative error",
    #         color=color_dict[method]
    #     ) 
    # else
    #     plot!(1:18,error_df[1:18,method],
    #         label=false,
    #         legend=:topleft,
    #         ylim=(10^-110,10^2),
    #         yscale=:log10,
    #         # title="relative error",
    #         color=color_dict[method],
    #         line=(2,:dash)
    #     ) 
    end
end
# for i in [1,5,6,7]
for i in [5,6,7]
    method = method_names[i]
    a,b = linear_fit(8:18,log10.(error_df[8:18,method]))
    plot!(8:35,10 .^(a .+ b.*[8:35...]),label=nothing,line=(1,:dash),yscale=:log10,color=color_dict[method])
end
xticks!([5,10,15,20,25,30,35], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}",L"2^{25}",L"2^{30}",L"2^{35}"])
plot!(legend = :bottom,legendfontsize = 6)
pb = plot(p3,p4,layout=(1,2),dpi=600) 
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/h_thread_fig2.png")




p5 = plot()
ylabel!("t (s)") 
for method in method_names[[1,7,6,4,5,2,3,8,9]]
    if occursin("thread",method) 
        if occursin("BF",method) || occursin("MF5",method); m=18; else; m=22; end; println(method," ",m)
        a,b = linear_fit(10:m,log10.(dt_df[10:m,method]))
        plot!(10:35,10 .^(a .+ b.*[10:35...]),label=nothing,line=(1,:dash),yscale=:log10,color=color_dict[method])
        println(1:m," ",dt_df[1:m,method])
        plot!(1:m,dt_df[1:m,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-4,10^8),
            yscale=:log10,
            # title="running time",
            color=color_dict[method]
        ) 
    end
end
xticks!([5,10,15,20,25,30,35], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}",L"2^{25}",L"2^{30}",L"2^{35}"])
plot!(legend = :bottomright,legendfontsize = 7)
p6 = plot()
ylabel!(L"$\delta$")
plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,linestyle=:dashdot,color="black") 
# plot!(1:18,pfplus_BF[1:18],label=L"p(F^+,F^-)",legend=:topleft,yscale=:log10,title="relative error",linestyle=:dashdot,color="black") 
a,b = linear_fit(8:18,log10.(pfplus_BF[8:18]))
plot!(8:35,10 .^(a .+ b.*[8:35...]),label=nothing,line=(1,:dash),yscale=:log10,color="black")
for i in [9,3,5,6,7,1]
    method = method_names[i]
    if occursin("thread",method) & !occursin("BF",method)
        hline!([eps(float_types[i])],color=color_dict[method],line=(0.3,:dashdot),label=false)
        annotate!(32.5, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
        plot!(1:18,pfplus_BF[1:18]./eps(float_types[i]),label=nothing,legend=:topleft,yscale=:log10,line=(0.3,:dashdot),color=color_dict[method]) 
        a,b = linear_fit(8:18,log10.(pfplus_BF[8:18]))
        plot!(8:35,10 .^(a .+ b.*[8:35...])./eps(float_types[i]),label=nothing,line=(0.3,:dashdot),yscale=:log10,color=color_dict[method])
        plot!(1:18,error_df[1:18,method],
            label=names_dict[method],
            legend=:topleft,
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
        ) 
    elseif method=="trick BF thread"
        hline!([eps(float_types[i])],color=color_dict[method],line=(0.3,:dashdot),label=false)
        annotate!(32.5, eps(float_types[i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
        setprecision(BigFloat,300) # because this is also done in quickscore algorithm
        plot!(1:18,pfplus_BF[1:18]./eps(float_types[i]),label=nothing,legend=:topleft,yscale=:log10,line=(0.3,:dashdot),color=color_dict[method]) 
        a,b = linear_fit(8:18,log10.(pfplus_BF[8:18]))
        plot!(8:35,10 .^(a .+ b.*[8:35...])./eps(float_types[i]),label=nothing,line=(0.3,:dashdot),yscale=:log10,color=color_dict[method])
        plot!(3:18,error_df[3:18,method],
            label=names_dict[method],
            legend=:topleft,
            ylim=(10^-100,10^10),
            yscale=:log10,
            # title="relative error",
            color=color_dict[method]
        ) 
    end
end
for i in [1,5,6,7]
    method = method_names[i]
    a,b = linear_fit(8:18,log10.(error_df[8:18,method]))
    plot!(8:35,10 .^(a .+ b.*[8:35...]),label=nothing,line=(1,:dash),yscale=:log10,color=color_dict[method])
end

xticks!([5,10,15,20,25,30,35], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}",L"2^{25}",L"2^{30}",L"2^{35}"])
plot!(legend = :bottom,legendfontsize = 6)
pc = plot(p5,p6,layout=(1,2),dpi=600) 
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/h_thread_fig3.png")

