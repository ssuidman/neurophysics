include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

file = XLSX.openxlsx("/Users/sam/Downloads/patient case.xlsx")
    data_alisa = file["Sheet1"][7:end,1:end]
close(file)

# 2867 true, 2502 true, 1731 true, 3568 false, 437 true, 3426 true, 3978 false, 1122 true, 1019 true,55 false, 3118 false, 987 false
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient404.mat"),"r")
    posterior_matlab = read(matfile, "pdiag")[:,1]
close(matfile)

sens = Float64.(data_alisa[2:end,3:end])[:,[1,2,3,5,6,8,9]]
sensneg = Float64.(data_alisa[2:end,3:end])[:,[4,7,10,11,12]]
previn = Float64.(data_alisa[2:end,2])
pfmin = 1 .- sens'
pfminneg = prod(1 .- sensneg', dims=1)[1,:]

m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
print("m=$m",'\n')
prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
prevend = vcat(1,previn) # prevend are the prevalences 1, P(d_1=1),... P(d_n=1) and is needed to multiply at the end
prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)


posteriors = DataFrame(nr = Int[1], Method = String["Matlab"], range = Vector{Float64}[[minimum(posterior_matlab),maximum(posterior_matlab)]], Posterior = Vector{BigFloat}[posterior_matlab], Pfplus = Vector{BigFloat}[[]])
methods = ["exp-sum-log","exp-sum-log BigFloat(all)","exp-sum-log prod1","prod","prod BigFloat(all)","prod prod1","prod BigFloat(pfmin)","exp-sum-log BigFloat(pfmin)"]
for (i,method) in enumerate(methods)
    println(method)
    pfplus, posterior, dt_array = quickscore(previn, pfmin, pfminneg, true, method)
    push!(posteriors, (nr = i+1,Method = method, range = [minimum(posterior),maximum(posterior)],Posterior = posterior, Pfplus = pfplus[:,1]))
end
pretty_table(posteriors[[1,2,3,9,5,6,8],["nr","Method","range","Posterior"]],alignment=:l)

x = log10.(abs.((Float64.((posteriors.Pfplus[8] - posteriors.Pfplus[9])[:,1]))))
y = log10.(abs.((Float64.((posteriors.Pfplus[3] - posteriors.Pfplus[6])[:,1]))))
println("Log(diff) posterior log/prod")
pretty_table(hcat(x,y),header=["BigFloat(pfmin)","BigFloat(all)"],alignment=:l)


