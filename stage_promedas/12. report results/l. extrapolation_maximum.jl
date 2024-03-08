include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

previn_dict,pfmin_dict,pfminneg_dict = prepare_patient_data("case 9"); m=size(pfmin_dict["1"],1); 
# previn,pfmin,pfminneg = prepare_patient_data("case 6"); m=size(pfmin); 
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn,pfmin[1:10,:],pfminneg,"trick BF thread");  
pfplus_matrix, pfplus1, P_joint, posterior1, dt = quickscore_report(previn,pfmin[1:10,:],pfminneg,"trick thread");  
pfplus_matrix, pfplus2, P_joint, posterior2, dt = quickscore_report(previn,pfmin[1:10,:],pfminneg,"trick MF2 thread");  

method_names = ["trick $(type)thread" for type in ["BF extreme ","MF6 ","MF5 ","MF4 ","MF3 ","MF2 ",""]]
float_types = [BigFloat,Float64x6,Float64x5,Float64x4,Float64x3,Float64x2,Float64]
names_dict = Dict("trick thread"=>"Float64","trick MF2 thread"=>"Float64x2","trick MF3 thread"=>"Float64x3","trick MF4 thread"=>"Float64x4","trick MF5 thread"=>"Float64x5","trick MF6 thread"=>"Float64x6","trick BF extreme thread"=>"BigFloat") 
eps_names_dict = Dict("trick BF thread"=>"BF","trick MF2"=>"64x2","trick MF2 thread"=>"64x2","trick MF3"=>"64x3","trick MF3 thread"=>"64x3","trick MF4 thread"=>"64x4","trick MF5 thread"=>"64x5","trick MF6 thread"=>"64x6","trick"=>"64","trick thread"=>"64") 
color_dict = Dict(method=>8-i for (i,method) in enumerate(method_names))
setprecision(BigFloat,350) # because this has also been done in quickscore algorithm 
eps_dict = Dict(method=>Float64(eps(float_types[i])) for (i,method) in enumerate(method_names))

error_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
dt_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
for k=[5]
    posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
    pfplus_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
    dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
    error_dict = Dict(method_name => zeros(m) for method_name in ["pfplus BF",method_names...])

    for i=1:m
        for method in method_names # important to first do 'BF', otherwise error_dict cannot be calculated
            println("k = $k")
            pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn_dict["$k"],pfmin_dict["$k"][1:i,:],pfminneg_dict["$k"],method)[[2,4,5]]; 
            error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF extreme thread"][i]) ./ posterior_dict["trick BF extreme thread"][i])
            if occursin("BF",method)
                error_dict["pfplus BF"][i] = pfplus[1,1]
            end 
        end
    end
    dt_df_dict["$k"] = DataFrame(dt_dict)
    error_df_dict["$k"] = DataFrame(error_dict)
end
# for k=1:length(previn_dict)
#     CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_dt.csv", dt_df_dict["$k"])    
#     CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_error.csv", error_df_dict["$k"])    
# end

error_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
dt_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
for k=1:length(previn_dict)
    dt_df_dict["$k"] = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_dt.csv", DataFrame)    
    error_df_dict["$k"] = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_error.csv", DataFrame)    
end

# Ik kan dit testen door 'h. thread.jl' te runnen t/m het moment dat dt_df en error_df gemaakt zijn 
extrapolation_df_dict = Dict{String,DataFrame}()
for k=1:length(previn_dict)
    error_df = error_df_dict["$k"]
    dt_df = dt_df_dict["$k"]
    extrapolation_df = DataFrame(
        "method"=>Vector{String}(undef,length(method_names)-1),
        "m_max"=>Vector{Int64}(undef,length(method_names)-1),
        "eps"=>Vector{Float64}(undef,length(method_names)-1),
        "pfplus"=>Vector{Float64}(undef,length(method_names)-1),
        "dt"=>Vector{Float64}(undef,length(method_names)-1),
        "error"=>Vector{Float64}(undef,length(method_names)-1),
        "error+1"=>Vector{Float64}(undef,length(method_names)-1)
        )
    for (i,type) in enumerate(["","MF2 ","MF3 ","MF4 ","MF5 ","MF6 "])
        method = "trick "*type*"thread"
        println("a: ",method)
        x = error_df[:,method]
        if nanmaximum(x) > 0.1
            nr = findall(x.==maximum(x[x.<0.1]))[1]
            error = x[nr]
            error_plus_1 = x[nr+1]
        else
            a_error,b_error = linear_fit(8:m,log10.(error_df[8:m,method])) 
            y_error = 10 .^(a_error .+ b_error.*[1:40...]) 
            nr = findall(y_error.==maximum(y_error[y_error.<0.1]))[1]
            error = y_error[nr]
            error_plus_1 = y_error[nr+1]
            # println(method," --->\t",nr,"\t",y_dt[nr])
        end
        println("b: ",nr)
        if nr<=m # check if the value <0.1 exists in the dt dataframe
            # println("1 ",method)
            dt = dt_df[nr,method]
            pfplus = error_df[nr,"pfplus BF"]
        else # if the value is not in the dt dataframe then do extrapolation 
            # println("3 ",method)
            a_dt,b_dt = linear_fit(8:m,log10.(dt_df[8:m,method])) 
            y_dt = 10 .^(a_dt .+ b_dt.*[1:40...]) 
            dt = y_dt[nr]
            a_pfplus,b_pfplus = linear_fit(8:m,log10.(error_df[8:m,"pfplus BF"])) 
            y_pfplus = 10 .^(a_pfplus .+ b_pfplus.*[1:40...]) 
            pfplus = y_pfplus[nr]
        end
        extrapolation_df[i,:] = [method,nr,eps_dict[method],pfplus,dt,error,error_plus_1]
        println(method," --->\t",nr,"\t",dt,"\t",error)
        extrapolation_df_dict["$k"] = extrapolation_df
    end 
end
# for k=1:length(previn_dict)
#     CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_extrapolation.csv", extrapolation_df_dict["$k"])    
# end

extrapolation_df_dict = Dict{String,DataFrame}()
for k=1:length(previn_dict)
    extrapolation_df_dict["$k"] = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_extrapolation.csv", DataFrame)    
end

gr() 
default(fontfamily="Times New Roman") 
# markers = [:circle, :utriangle, :dtriangle, :square, :diamond, :pentagon, :hexagon, :star4, :star5, :star6]
# markers = [:circle,:utriangle,:utriangle,:utriangle,:circle,:utriangle,:utriangle,:circle,:utriangle,:utriangle]
markers = [:circle for i=1:10]
p1 = plot(21:32,fill(10^-4,12),  fillrange=fill(10^6.1,12),fillalpha=0.15,color="red",label=false,yscale=:log10,xlims=(0,31),ylims=(10^-4,10^6),xlabel="max #iterations for "*latexstring("\\delta<0.1"),ylabel=L"t(s) \cdot \frac{1000}{n}",xticks=([5:5:30...],[latexstring("2^{$k}") for k in 5:5:30]),yticks=10.0.^[-4:6...],legend=:topleft)
p2 = plot(21:32,fill(10^-4,12),  fillrange=fill(10^6.1,12),fillalpha=0.15,color="red",label=false,yscale=:log10,xlims=(0,31),ylims=(10^-4,10^6),xlabel="max #iterations for "*latexstring("\\delta<0.1"),ylabel=L"t(s) \cdot \frac{1000}{n}",xticks=([5:5:30...],[latexstring("2^{$k}") for k in 5:5:30]),yticks=10.0.^[-4:6...],legend=:topleft)
p3 = plot(21:32,fill(10^-1,12),  fillrange=fill(10^15,12), fillalpha=0.15,color="red",label=false,yscale=:log10,xlims=(0,31),ylims=(10^0,10^13),xlabel="max #iterations for "*latexstring("\\delta<0.1"),ylabel=L"\tau"                     ,xticks=([5:5:35...],[latexstring("2^{$k}") for k in 5:5:30]),yticks=10.0.^[0:13...],legend=:topleft)
p4 = plot(21:32,fill(10^-100,12),fillrange=fill(10^15,12), fillalpha=0.15,color="red",label=false,yscale=:log10,xlims=(0,31),ylims=(10^-100,10^0),xlabel="max #iterations for "*latexstring("\\delta<0.1"),ylabel=L"p(F^+,F^-)"             ,xticks=([5:5:35...],[latexstring("2^{$k}") for k in 5:5:30]),yticks=10.0.^[-100:10:0...],legend=:bottomleft)
for k in 1:length(previn_dict)
    n = length(previn_dict["$k"])
    extrapolation_df = extrapolation_df_dict["$k"]
    for i in 6:-1:1
        method = extrapolation_df[i,"method"]
        if k==1
            scatter!(p1,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]],label=names_dict[method],marker=markers[k],color=color_dict[method])
            scatter!(p2,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]/n*1000],label=names_dict[method],marker=markers[k],color=color_dict[method])
            scatter!(p3,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"pfplus"]]/[eps_dict[method]],label=names_dict[method],marker=markers[k],color=color_dict[method])
            scatter!(p4,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"pfplus"]],label=names_dict[method],marker=markers[k],color=color_dict[method])
            hline!(p4,[eps(float_types[8-i])],color=color_dict[method],line=(0.3,:dashdot),label=false)
            annotate!(p4,28.5, eps(float_types[8-i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
        else
            scatter!(p1,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]],label=false,marker=markers[k],color=color_dict[method])
            scatter!(p2,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]/n*1000],label=false,marker=markers[k],color=color_dict[method])
            scatter!(p3,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"pfplus"]]/[eps_dict[method]],label=false,marker=markers[k],color=color_dict[method])
            scatter!(p4,[extrapolation_df[i,"m_max"]],[extrapolation_df[i,"pfplus"]],label=false,marker=markers[k],color=color_dict[method])
            hline!(p4,[eps(float_types[8-i])],color=color_dict[method],line=(0.3,:dashdot),label=false)
            annotate!(p4,28.5, eps(float_types[8-i]), text(latexstring("\\textbf{\\epsilon_{\\textbf{$(eps_names_dict[method])}}}"), :bottom, 7, :black))
        end
    end
end
vline!(p1,[m+1],color="red",label=false)
annotate!(p1,16.75,10^-2,text("measured",7,"black"))
annotate!(p1,26.25,10^-2,text("extrapolated",7,"red"))
vline!(p2,[m+1],color="red",label=false)
annotate!(p2,16.75,10^-2,text("measured",7,"black"))
annotate!(p2,26.25,10^-2,text("extrapolated",7,"red"))
hline!(p2,[60],color="black",line=(0.3,:dash),label=false)
annotate!(p2,15, [60], text("minute", :bottom, 7, :black))
hline!(p2,[60*60],color="black",line=(0.3,:dash),label=false)
annotate!(p2,15, [60*60], text("hour", :bottom, 7, :black))
hline!(p2,[60*60*24],color="black",line=(0.3,:dash),label=false)
annotate!(p2,15, [60*60*24], text("day", :bottom, 7, :black))
vline!(p3,[m+1],color="red",label=false)
annotate!(p3,16.75,10^2,text("measured",7,"black"))
annotate!(p3,26.25,10^2,text("extrapolated",7,"red"))
vline!(p4,[m+1],color="red",label=false)
annotate!(p4,16.75,10^-8,text("measured",7,"black"))
annotate!(p4,26.25,10^-8,text("extrapolated",7,"red"))
# p = plot(p1,p2,layout=(1,2),dpi=600)
p = plot(p2,p3,layout=(1,2),dpi=600)
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/l_extrapolation_maximum_fig1.png")
p = plot(p2,p4,layout=(1,2),dpi=600)
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/l_extrapolation_maximum_fig2.png")

# Calculate on avarage how much longer e.g. MF6 is than MF5
all_running_times = reduce(vcat,[Matrix(dt_df_dict["$k"]) for k=1:10]);
MF6_divided_MF5 = mean(all_running_times[:,6]./all_running_times[:,5])
MF5_divided_MF4 = mean(all_running_times[:,5]./all_running_times[:,4])
MF4_divided_MF3 = mean(all_running_times[:,4]./all_running_times[:,3])
MF3_divided_MF2 = mean(all_running_times[:,3]./all_running_times[:,2])
MF2_divided_64 = mean(all_running_times[:,2]./all_running_times[:,7])


