include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# 2867 true, 2502 true, 1731 true, 3568 false, 437 true, 3426 true, 3978 false, 1122 true, 1019 true,55 false, 3118 false, 987 false
matrix_alisa = Dict{String,Matrix{Any}}()
file = XLSX.openxlsx("/Users/sam/Downloads/patient case1.xlsx")
    matrix_alisa["case 1"] = file["Sheet1"][1:end, 1:end]
    matrix_alisa["case 2"] = file["Sheet2"][1:end, 1:end]
    matrix_alisa["case 3"] = file["Sheet3"][1:end, 1:end]
close(file)

matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient404.mat"),"r")
    posterior_matlab = read(matfile, "pdiag")[:,1]
    P_joint_matlab = read(matfile, "P_joint")[:,1]
    dt_matlab = read(matfile,"dt")
close(matfile)




############ HIER ONDER EN BOVEN VERDER MET DICTIONARIES EN CASES ZODAT IK SHEETS GOED KAN INLADEN ###############
matrix_alisa["case 1"][4,1:14]
matrix_alisa["case 3"][4,1:15]
data_alisa = Dict(
    "case 1"=>DataFrame(Float64.(matrix_alisa["case 1"][8:end,:]),matrix_alisa["case 1"][6,:]),
    "case 3"=>DataFrame(Float64.(matrix_alisa["case 3"][8:end,:]),matrix_alisa["case 3"][6,:])
)
Matrix(data_alisa[case][:,idx_sens])[:,[1,2,3,5,6,8,9]]

diagnosis_nr, P_joint_hp_alisa, P_joint_alisa, posterior_hp_alisa, posterior_alisa = Dict{String,Vector{Int64}}(), Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}()
for case in keys(data_alisa)
    diagnosis_nr[case] = Int64.(data_alisa[case][:,"diagnosis ID number"])
    P_joint_hp_alisa[case] = data_alisa[case][:,"Pjoint (high precision)"]
    P_joint_alisa[case] = data_alisa[case][:,"Pjoint(standart double precision)"]
    posterior_hp_alisa[case] = data_alisa[case][:,"Pdiag (high precision)"]
    posterior_alisa[case] = data_alisa[case][:,"Pdiag(standart double precision)"]

    idx_sens = names(data_alisa[case])[startswith.(names(data_alisa[case]),"Sensitivity")]
    sens = Matrix(data_alisa[case][:,idx_sens])[:,[1,2,3,5,6,8,9]]
    sensneg = Matrix(data_alisa[case][:,idx_sens])[:,[4,7,10,11,12]]
    previn = data_alisa[case][:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]
end




posteriors = DataFrame(
    nr = Int[1,2,3], 
    Method = String["Matlab","Alisa normal","Alisa hp"], 
    Posterior_range = Vector{Float64}[[minimum(posterior_matlab),maximum(posterior_matlab)],[minimum(posterior_alisa),maximum(posterior_alisa)],[minimum(posterior_hp_alisa),maximum(posterior_hp_alisa)]], 
    Posterior = Vector{BigFloat}[posterior_matlab,posterior_alisa,posterior_hp_alisa], 
    P_joint_range = Vector{Float64}[[minimum(P_joint_matlab),maximum(P_joint_matlab)],[minimum(P_joint_alisa),maximum(P_joint_alisa)],[minimum(P_joint_hp_alisa),maximum(P_joint_hp_alisa)]], 
    P_joint = Vector{BigFloat}[P_joint_matlab,P_joint_alisa,P_joint_hp_alisa],
    time = Float64[dt_matlab,0,0]
)
methods = ["exp-sum-log","exp-sum-log BF(all)","prod","prod BF(all)","prod BF(pfmin)"]#,"exp-sum-log BF(pfmin)"]
for (i,method) in enumerate(methods)
    println(method)
    pfplus, posterior, dt = quickscore(previn, pfmin, pfminneg, true, method)
    push!(posteriors, (nr = i+3,
        Method = method, 
        Posterior_range = [minimum(posterior),maximum(posterior)],
        Posterior = posterior, 
        P_joint_range = [minimum(pfplus[2:end,1].*previn),maximum(pfplus[2:end,1].*previn)], 
        P_joint = pfplus[2:end,1].*previn,
        time = dt))
end
pretty_table(posteriors[:,["nr","Method","P_joint_range","Posterior_range","time"]],alignment=:l)


diff_posterior, diff_P_joint = Vector{Float64}[], Vector{Float64}[]
for i=eachindex(posteriors.Method)
    push!(diff_posterior,Float64.(maximum(abs,reduce(hcat,posteriors.Posterior) .- posteriors.Posterior[i],dims=1)[1,:]))
    push!(diff_P_joint,Float64.(maximum(abs,reduce(hcat,posteriors.P_joint) .- posteriors.P_joint[i],dims=1)[1,:]))
end
idx = 1:8
println("############################################## max-abs diffference of P_joint between methods ##############################################")
pretty_table(hcat(reduce(hcat,diff_P_joint)[idx,idx],posteriors.time),row_names=posteriors.Method[idx],header=vcat(posteriors.Method[idx],"time"),alignment=:l,formatters=ft_printf("%1.0e"))
println("########################################### max-abs diffference of Posteriors between methods ###########################################")
pretty_table(hcat(reduce(hcat,diff_posterior)[idx,idx],posteriors.time),row_names=posteriors.Method[idx],header=vcat(posteriors.Method[idx],"time"),alignment=:l)


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





