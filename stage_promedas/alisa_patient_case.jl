include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

matrix_alisa = Dict{String,Matrix{Any}}()
file = XLSX.openxlsx("/Users/sam/Downloads/patient case1.xlsx")
    matrix_alisa["case 1"] = file["Sheet1"][1:end, 1:end]
    matrix_alisa["case 2"] = file["Sheet2"][1:end, 1:end]
    matrix_alisa["case 3"] = file["Sheet3"][1:end, 1:end]
close(file)

matrix_alisa["case 2"][7,:]
matrix_alisa["case 3"][7,:]

posterior_matlab, P_joint_matlab, dt_matlab = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient404.mat"),"r")
    P_joint_matlab["case 1"] = read(matfile, "P_joint")[:,1]
    posterior_matlab["case 1"] = read(matfile, "pdiag")[:,1]
    dt_matlab["case 1"] = read(matfile,"dt")
    P_joint_matlab["case 2"] = read(matfile, "P_joint")[:,1]
    posterior_matlab["case 2"] = read(matfile, "pdiag")[:,1]
    dt_matlab["case 2"] = read(matfile,"dt")
    P_joint_matlab["case 3"] = read(matfile, "P_joint")[:,1]
    posterior_matlab["case 3"] = read(matfile, "pdiag")[:,1]
    dt_matlab["case 3"] = read(matfile,"dt")
    # for case_nr=1:3
    #     posterior_matlab["case $case_nr"] = read(matfile, "pdiag_case_$case_nr")[:,1]
    #     P_joint_matlab["case $case_nr"] = read(matfile, "P_joint_case_$case_nr")[:,1]
    #     dt_matlab["case $case_nr"] = read(matfile,"dt_case_$case_nr")
    # end
close(matfile)

function posteriors_func(case,methods,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)
    data_alisa = DataFrame(Float64.(matrix_alisa[case][8:end,:]),matrix_alisa[case][6,:])
    # diagnosis_nr = Int64.(data_alisa[:,"diagnosis ID number"])
    P_joint_hp_alisa = data_alisa[:,"Pjoint (high precision)"]
    P_joint_alisa = data_alisa[:,"Pjoint(standart double precision)"]
    posterior_hp_alisa = data_alisa[:,"Pdiag (high precision)"]
    posterior_alisa = data_alisa[:,"Pdiag(standart double precision)"]

    idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
    all_sens = Matrix(data_alisa[:,idx_sens])
    pos_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'true'"
    neg_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'false'"
    sens = all_sens[:,pos_idx]
    sensneg = all_sens[:,neg_idx]
    previn = data_alisa[:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]
    
    posteriors = DataFrame(
        nr = Int[1,2,3], 
        Method = String["Matlab","Alisa normal","Alisa hp"], 
        # Posterior_range = Vector{Float64}[[minimum(posterior_matlab[case]),maximum(posterior_matlab[case])],[minimum(posterior_alisa),maximum(posterior_alisa)],[minimum(posterior_hp_alisa),maximum(posterior_hp_alisa)]], 
        Posterior_min = Float64[minimum(posterior_matlab[case]),minimum(posterior_alisa),minimum(posterior_hp_alisa)], 
        Posterior_max = Float64[maximum(posterior_matlab[case]),maximum(posterior_alisa),maximum(posterior_hp_alisa)], 
        Posterior = Vector{BigFloat}[posterior_matlab[case],posterior_alisa,posterior_hp_alisa], 
        # P_joint_range = Vector{Float64}[[minimum(P_joint_matlab[case]),maximum(P_joint_matlab[case])],[minimum(P_joint_alisa),maximum(P_joint_alisa)],[minimum(P_joint_hp_alisa),maximum(P_joint_hp_alisa)]], 
        P_joint_min = Float64[minimum(P_joint_matlab[case]),minimum(P_joint_alisa),minimum(P_joint_hp_alisa)], 
        P_joint_max = Float64[maximum(P_joint_matlab[case]),maximum(P_joint_alisa),maximum(P_joint_hp_alisa)], 
        P_joint = Vector{BigFloat}[P_joint_matlab[case],P_joint_alisa,P_joint_hp_alisa],
        time = Float64[dt_matlab[case],0,0]
    )

    for (i,method) in enumerate(methods)
        println(method)
        pfplus, posterior, dt = quickscore(previn, pfmin, pfminneg, true, method)
        push!(posteriors, (nr = i+3,
            Method = method, 
            # Posterior_range = [minimum(posterior),maximum(posterior)],
            Posterior_min = minimum(posterior),
            Posterior_max = maximum(posterior),
            Posterior = posterior, 
            # P_joint_range = [minimum(pfplus[2:end,1].*previn),maximum(pfplus[2:end,1].*previn)], 
            P_joint_min = minimum(pfplus[2:end,1].*previn),
            P_joint_max = maximum(pfplus[2:end,1].*previn),
            P_joint = pfplus[2:end,1].*previn,
            time = dt))
    end

    return posteriors
end

methods_1 = ["exp-sum-log","exp-sum-log BF(all)","prod","prod BF(all)","prod BF(pfmin)"]
methods_2 = ["exp-sum-log","exp-sum-log BF(all)","prod","prod BF(all)","prod BF(pfmin)","exp-sum-log BF(pfmin)"]
posteriors = Dict{String,DataFrame}()
posteriors["case 1"] = posteriors_func("case 1",methods_1,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)
posteriors["case 2"] = posteriors_func("case 2",methods_1,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)
posteriors["case 3"] = posteriors_func("case 3",methods_1,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)

matrix_alisa["case 1"][7,:]
matrix_alisa["case 2"][7,:]
matrix_alisa["case 3"][7,:]



function print_posteriors(posteriors,lim1,lim2)    
    diff_posterior, diff_P_joint = Vector{Float64}[], Vector{Float64}[]
    for i=eachindex(posteriors.Method)
        push!(diff_posterior,Float64.(maximum(abs,reduce(hcat,posteriors.Posterior) .- posteriors.Posterior[i],dims=1)[1,:]))
        push!(diff_P_joint,Float64.(maximum(abs,reduce(hcat,posteriors.P_joint) .- posteriors.P_joint[i],dims=1)[1,:]))
    end
    hl1 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] >= 0 && d[i, j] <= 1) : false, crayon"green")
    hl2 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] <= 0 || d[i, j] >= 1) : false, crayon"red")
    hl3 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim1 && d[i,j]!=0) : false, crayon"green")
    hl4 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim2 && d[i,j]!=0) : false, crayon"green")

    println("############################################## range of P_joint and Posterior between methods ##############################################")
    pretty_table(posteriors[:,["nr","Method","P_joint_min","P_joint_max","Posterior_min","Posterior_max"]],alignment=:l,highlighters=(hl1,hl2))
    println("############################################## max-abs diffference of P_joint between methods ##############################################")
    pretty_table(hcat(reduce(hcat,diff_P_joint),posteriors.time),row_names=posteriors.Method,header=vcat(posteriors.Method,"time"),alignment=:l,highlighters=hl3) # ,formatters=ft_printf("%1.0e")
    println("########################################### max-abs diffference of Posterior between methods ###########################################")
    pretty_table(hcat(reduce(hcat,diff_posterior),posteriors.time),row_names=posteriors.Method,header=vcat(posteriors.Method,"time"),alignment=:l,highlighters=hl4)
end 

println("############################################## case 1 (m=7) ##############################################")
print_posteriors(posteriors["case 1"],1e-17,1e-6)
println("############################################## case 2 (m=8) ##############################################")
print_posteriors(posteriors["case 2"][2:end,:],1e-17,1e-6)
println("############################################## case 3 (m=9) ##############################################")
print_posteriors(posteriors["case 3"][2:end,:],1e-17,1e-6)



# println("normal BigFloat")
# @btime prod(BigFloat.(pfmin),dims=1);
# println("normal BigFloat with specification of precision = 256 bits")
# @btime prod(BigFloat.(pfmin,precision=256),dims=1);
# println("precision of 100 bits")
# @btime prod(BigFloat.(pfmin,precision=100),dims=1);
# println("same precision = 53 bits as Float64")
# @btime prod(BigFloat.(pfmin,precision=53),dims=1);
# println("normal Float64 (no BigFloat at all)")
# @btime prod(pfmin,dims=1);





