include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

# previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin);
previn_dict,pfmin_dict,pfminneg_dict = prepare_patient_data("case 9"); 

method_names = ["trick $(type)thread" for type in ["","MF2 ","MF3 ","MF4 ","MF5 ","BF "]]
float_types = [Float64,Float64x2,Float64x3,Float64x4,Float64x5,BigFloat]
names_dict = Dict("trick thread"=>"Float64","trick MF2 thread"=>"Float64x2","trick MF3 thread"=>"Float64x3","trick MF4 thread"=>"Float64x4","trick MF5 thread"=>"Float64x5","trick BF thread"=>"BigFloat") 
color_dict = Dict("trick thread"=>1,"trick MF2 thread"=>2,"trick MF3 thread"=>3,"trick MF4 thread"=>4,"trick MF5 thread"=>5,"trick BF thread"=>6) 

posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
error_dict = Dict(method_name => zeros(m) for method_name in method_names)
pfplus_BF = zeros(BigFloat,m)
# For the BF methods only up to m=18
for i=1:m
    for method in method_names
        # if occursin("BF",method) && i>18
        #     posterior_dict[method][i],dt_dict[method][i] = zeros(n), 0.0 
        # elseif occursin("BF",method) && i<=18
        #     pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[2,4,5]]; 
        #     error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
        #     pfplus_BF[i] = pfplus[1,1]
        # elseif !occursin("BF",method) && i<=18
        #     posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
        #     error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
        # else
        # posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[4,5]]; 
        # end
        pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn,pfmin[1:i,:],pfminneg,method)[[2,4,5]]; 
        error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
        pfplus_BF[i] = pfplus[1,1]
    end
end
dt_df = DataFrame(dt_dict)[:,method_names[1:end]]
error_df = DataFrame(error_dict)[:,method_names[1:end]]
# Then calculate BigFloat pfplus 
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_time.csv",DataFrame(dt_dict)[:,["trick MF5 thread"]])
# CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_error.csv", DataFrame(error_dict)[:,["trick MF5 thread"]])
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_pfplus_BF.csv",pfplus_BF)

plot(1:m,error_df[1:m,method_names[5]],yscale=:log10)

dt_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_time.csv",DataFrame)
error_df = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_error.csv",DataFrame)
pfplus_BF = readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/thread_pfplus_BF.csv")



# Ik kan dit testen door 'h. thread.jl' te runnen t/m het moment dat dt_df en error_df gemaakt zijn 
setprecision(BigFloat,300) # because this has also been done in quickscore algorithm 
eps_dict = Dict("trick thread"=>eps(Float64),"trick MF2 thread"=>eps(Float64x2),"trick MF3 thread"=>eps(Float64x3),"trick MF4 thread"=>eps(Float64x4),"trick MF5 thread"=>eps(Float64x5),"trick BF thread"=>eps(BigFloat))
extrapolation_df = DataFrame("method"=>Vector{String}(undef,6),"eps"=>Vector{Float64}(undef,6),"m_max"=>Vector{Int64}(undef,6),"dt"=>Vector{Float64}(undef,6),"error"=>Vector{Float64}(undef,6),"error+1"=>Vector{Float64}(undef,6))
for (i,type) in enumerate(["","MF2 ","MF3 ","MF4 ","MF5 ","BF "])
    method = "trick "*type*"thread"
    x = error_df[1:18,method]
    if maximum(x) > 0.1
        nr = findall(x.==maximum(x[x.<0.1]))[1]
        error = x[nr]
        error_plus_1 = x[nr+1]
    else
        a_error,b_error = linear_fit(8:18,log10.(error_df[8:18,method])) 
        y_error = 10 .^(a_error .+ b_error.*[1:35...]) 
        nr = findall(y_error.==maximum(y_error[y_error.<0.1]))[1]
        error = y_error[nr]
        error_plus_1 = y_error[nr+1]
        # println(method," --->\t",nr,"\t",y_dt[nr])
    end
    if nr<=22 # check if the value <0.1 exists in the dt dataframe
        if dt_df[nr,method]!=0 # check if the value is non-zero 
            println("1 ",method)
            dt = dt_df[nr,method]
        else # if the value is zero then do extrapolation
            println("2 ",method)
            a_dt,b_dt = linear_fit(8:18,log10.(dt_df[8:18,method])) 
            y_dt = 10 .^(a_dt .+ b_dt.*[1:35...]) 
            dt = y_dt[nr]
        end
    else # if the value is not in the dt dataframe then do extrapolation 
        println("3 ",method)
        a_dt,b_dt = linear_fit(8:18,log10.(dt_df[8:18,method])) 
        y_dt = 10 .^(a_dt .+ b_dt.*[1:35...]) 
        dt = y_dt[nr]
    end
    extrapolation_df[i,:] = [method,eps_dict[method],nr,dt,error,error_plus_1]
    println(method," --->\t",nr,"\t",dt,"\t",error)
end 
extrapolation_df

p = plot()
for i in 6:-1:1
    scatter!([extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]],yscale=:log10,label=names_dict[method_names[i]],legend=:topleft)
end
display(p)


