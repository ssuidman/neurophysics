include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")
include("0. prepare_patient_case_9.jl")

# previn, pfmin, pfminneg = prepare_patient_data("case 5"); (m,n)=size(pfmin);
previn_dict,pfmin_dict,pfminneg_dict = prepare_patient_data("case 9"); m=size(pfmin_dict["1"],1); 

method_names = ["trick $(type)thread" for type in ["BF ","MF5 ","MF4 ","MF3 ","MF2 ",""]]
float_types = [BigFloat,Float64x5,Float64x4,Float64x3,Float64x2,Float64]
names_dict = Dict("trick thread"=>"Float64","trick MF2 thread"=>"Float64x2","trick MF3 thread"=>"Float64x3","trick MF4 thread"=>"Float64x4","trick MF5 thread"=>"Float64x5","trick BF thread"=>"BigFloat") 
color_dict = Dict("trick thread"=>1,"trick MF2 thread"=>2,"trick MF3 thread"=>3,"trick MF4 thread"=>4,"trick MF5 thread"=>5,"trick BF thread"=>6) 
setprecision(BigFloat,300) # because this has also been done in quickscore algorithm 
eps_dict = Dict("trick thread"=>eps(Float64),"trick MF2 thread"=>eps(Float64x2),"trick MF3 thread"=>eps(Float64x3),"trick MF4 thread"=>eps(Float64x4),"trick MF5 thread"=>eps(Float64x5),"trick BF thread"=>eps(BigFloat))

error_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
dt_df_dict = Dict("$k"=>DataFrame() for k=1:length(previn_dict))
for k=1:length(previn_dict)
    posterior_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
    pfplus_dict = Dict(method_name => [Vector{float_types[i]}() for j in 1:m] for (i,method_name) in enumerate(method_names))
    dt_dict = Dict(method_name => zeros(m) for method_name in method_names)
    error_dict = Dict(method_name => zeros(m) for method_name in ["pfplus BF",method_names...])

    for i=1:m
        for method in method_names # important to first do 'BF', otherwise error_dict cannot be calculated
            println("k = $k")
            pfplus,posterior_dict[method][i],dt_dict[method][i] = quickscore_report(previn_dict["$k"],pfmin_dict["$k"][1:i,:],pfminneg_dict["$k"],method)[[2,4,5]]; 
            error_dict[method][i] = maximum(abs,(posterior_dict[method][i] .- posterior_dict["trick BF thread"][i]) ./ posterior_dict["trick BF thread"][i])
            if occursin("BF",method)
                error_dict["pfplus BF"][i] = pfplus[1,1]
            end 
        end
    end
    dt_df_dict["$k"] = DataFrame(dt_dict)
    error_df_dict["$k"] = DataFrame(error_dict)
end
# for k=1:length(previn_dict)
    # CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_dt.csv", dt_df_dict["$k"])    
    # CSV.write("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_error.csv", error_df_dict["$k"])    
# end

error_df_dict = Dict("$k"=>DataFrame() for k=1:10)
dt_df_dict = Dict("$k"=>DataFrame() for k=1:10)
for k=1:length(previn_dict)
    dt_df_dict["$k"] = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_dt.csv", DataFrame)    
    error_df_dict["$k"] = CSV.read("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/extrapolation_random_patients/patient_$(k)_error.csv", DataFrame)    
end


# Ik kan dit testen door 'h. thread.jl' te runnen t/m het moment dat dt_df en error_df gemaakt zijn 
error_df = error_df_dict["4"]
dt_df = dt_df_dict["4"]
extrapolation_df = DataFrame("method"=>Vector{String}(undef,5),"m_max"=>Vector{Int64}(undef,5),"eps"=>Vector{Float64}(undef,5),"pfplus"=>Vector{Float64}(undef,5),"dt"=>Vector{Float64}(undef,5),"error"=>Vector{Float64}(undef,5),"error+1"=>Vector{Float64}(undef,5))
for (i,type) in enumerate(["","MF2 ","MF3 ","MF4 ","MF5 "])
    method = "trick "*type*"thread"
    println("a: ",method)
    x = error_df[:,method]
    if nanmaximum(x) > 0.1
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
    println("b: ",nr)
    if nr<=m # check if the value <0.1 exists in the dt dataframe
        # println("1 ",method)
        dt = dt_df[nr,method]
        pfplus = error_df[nr,"pfplus BF"]
    else # if the value is not in the dt dataframe then do extrapolation 
        # println("3 ",method)
        a_dt,b_dt = linear_fit(8:18,log10.(dt_df[8:18,method])) 
        y_dt = 10 .^(a_dt .+ b_dt.*[1:35...]) 
        dt = y_dt[nr]
        a_pfplus,b_pfplus = linear_fit(8:18,log10.(error_df[8:18,"pfplus BF"])) 
        y_pfplus = 10 .^(a_pfplus .+ b_pfplus.*[1:35...]) 
        pfplus = y_pfplus[nr]
    end
    extrapolation_df[i,:] = [method,nr,eps_dict[method],pfplus,dt,error,error_plus_1]
    println(method," --->\t",nr,"\t",dt,"\t",error)
end 
extrapolation_df

p = plot()
for i in 6:-1:1
    scatter!([extrapolation_df[i,"m_max"]],[extrapolation_df[i,"dt"]],yscale=:log10,label=names_dict[method_names[i]],legend=:topleft)
end
display(p)


