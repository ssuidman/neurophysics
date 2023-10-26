include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# 2867 true, 2502 true, 1731 true, 3568 false, 437 true, 3426 true, 3978 false, 1122 true, 1019 true,55 false, 3118 false, 987 false
file = XLSX.openxlsx("/Users/sam/Downloads/patient case.xlsx")
    data_alisa = file["Sheet1"][7:end,1:end]
close(file)
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient404.mat"),"r")
    posterior_matlab = read(matfile, "pdiag")[:,1]
    P_joint_matlab = read(matfile, "P_joint")[:,1]
close(matfile)

P_joint_hp_alisa = Float64.(data_alisa[2:end,3])
P_joint_alisa = Float64.(data_alisa[2:end,4])
posterior_hp_alisa = Float64.(data_alisa[2:end,5])
posterior_alisa = Float64.(data_alisa[2:end,6])

sens = Float64.(data_alisa[2:end,7:end])[:,[1,2,3,5,6,8,9]]
sensneg = Float64.(data_alisa[2:end,7:end])[:,[4,7,10,11,12]]
previn = Float64.(data_alisa[2:end,2])
pfmin = 1 .- sens'
pfminneg = prod(1 .- sensneg', dims=1)[1,:]

posteriors = DataFrame(nr = Int[1,2,3], Method = String["Matlab","Alisa normal","Alisa hp"], Posterior_range = Vector{Float64}[[minimum(posterior_matlab),maximum(posterior_matlab)],[minimum(posterior_alisa),maximum(posterior_alisa)],[minimum(posterior_hp_alisa),maximum(posterior_hp_alisa)]], Posterior = Vector{BigFloat}[posterior_matlab,posterior_alisa,posterior_hp_alisa], P_joint_range = Vector{Float64}[[minimum(P_joint_matlab),maximum(P_joint_matlab)],[minimum(P_joint_alisa),maximum(P_joint_alisa)],[minimum(P_joint_hp_alisa),maximum(P_joint_hp_alisa)]], P_joint = Vector{BigFloat}[P_joint_matlab,P_joint_alisa,P_joint_hp_alisa])
methods = ["exp-sum-log","exp-sum-log BF(all)","prod","prod BF(all)","prod BF(pfmin)"]#,"exp-sum-log BF(pfmin)"]
for (i,method) in enumerate(methods)
    println(method)
    pfplus, posterior, dt_array = quickscore(previn, pfmin, pfminneg, true, method)
    push!(posteriors, (nr = i+3,Method = method, Posterior_range = [minimum(posterior),maximum(posterior)],Posterior = posterior, P_joint_range = [minimum(pfplus[2:end,1].*previn),maximum(pfplus[2:end,1].*previn)], P_joint = pfplus[2:end,1].*previn))
end
pretty_table(posteriors[:,["nr","Method","Posterior_range","Posterior"]],alignment=:l)
pretty_table(posteriors[:,["nr","Method","P_joint_range","P_joint"]],alignment=:l)
pretty_table(posteriors[:,["nr","Method","P_joint_range","Posterior_range"]],alignment=:l)

diff_posterior, diff_P_joint = Vector{Float64}[], Vector{Float64}[]
for i=eachindex(posteriors.Method)
    push!(diff_posterior,Float64.(maximum(abs,reduce(hcat,posteriors.Posterior) .- posteriors.Posterior[i],dims=1)[1,:]))
    push!(diff_P_joint,Float64.(maximum(abs,reduce(hcat,posteriors.P_joint) .- posteriors.P_joint[i],dims=1)[1,:]))
end
idx = 1:8
posteriors
println("############################################## max-abs diffference of P_joint between methods ##############################################")
pretty_table(reduce(hcat,diff_P_joint)[idx,idx],row_names=posteriors.Method[idx],header=posteriors.Method[idx],alignment=:l)
println("########################################### max-abs diffference of Posteriors between methods ###########################################")
pretty_table(reduce(hcat,diff_posterior)[idx,idx],row_names=posteriors.Method[idx],header=posteriors.Method[idx],alignment=:l)


println("normal BigFloat")
@btime prod(BigFloat.(pfmin),dims=1);
println("normal BigFloat with specification of precision = 256 bits")
@btime prod(BigFloat.(pfmin,precision=256),dims=1);
println("precision of 100 bits")
@btime prod(BigFloat.(pfmin,precision=100),dims=1);
println("same precision = 53 bits as Float64")
@btime prod(BigFloat.(pfmin,precision=53),dims=1);
println("normal Float64 (no BigFloat at all)")
@btime prod(pfmin,dims=1);





