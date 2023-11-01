include("../packages.jl")
include("../quickscore_algorithm.jl")

# Data preparation
matrix_alisa = Dict{String,Matrix{Any}}()
file = XLSX.openxlsx("/Users/sam/Downloads/patient cases.xlsx")
    matrix_alisa["case 1"] = file["Sheet1"][1:end, 1:end]
close(file)
case = "case 1"
data_alisa = DataFrame(Float64.(matrix_alisa[case][8:end,:]),matrix_alisa[case][6,:])
idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
all_sens = Matrix(data_alisa[:,idx_sens])
pos_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'true'"
neg_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'false'"
sens = all_sens[:,pos_idx]
sensneg = all_sens[:,neg_idx]
previn = data_alisa[:,"prevalence"]
pfmin = 1 .- sens'
pfminneg = prod(1 .- sensneg', dims=1)[1,:]
pfplus, posterior, dt = quickscore(previn, pfmin, pfminneg, true, "prod BF(pfmin)")

# Write important variables to csv-files 
writedlm("cplusplus_test/patient404_case_1_previn.csv", previn, ',')
writedlm("cplusplus_test/patient404_case_1_pfmin.csv", pfmin, ',')
writedlm("cplusplus_test/patient404_case_1_pfminneg.csv", pfminneg, ',')
writedlm("cplusplus_test/patient404_case_1_posterior.csv", posterior, ',')

# Retrieve variables from csv file 
previn = readdlm("cplusplus_test/patient404_case_1_previn.csv",',')[:,1]
pfmin = readdlm("cplusplus_test/patient404_case_1_pfmin.csv",',')
pfminneg = readdlm("cplusplus_test/patient404_case_1_pfminneg.csv",',')[:,1]
# posterior = readdlm("cplusplus_test/patient404_case_1_posterior.csv",',',BigFloat)[:,1]
posterior_ = readdlm("cplusplus_test/patient404_case_1_posterior.csv",',',BigFloat)[:,1]

m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
prevend = vcat(1,previn) # prevend are the prevalences 1, P(d_1=1),... P(d_n=1) and is needed to multiply at the end
prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
pfplus, pfplus_, myset = zeros(Float64,n+1), zeros(Float64,n+1), Vector{Int64}()
t = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements
        pfplus = pfplus .+ ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset,:],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
        # pfplus = pfplus .+ ((-1)^length(myset)) .* prod(1e-50 .+ (prod(BigFloat.(pfmin[myset, :]),dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
        # for j=1:n+1
        #     sumExp = 0.0
        #     for k=1:n
        #         product = 1.0
        #         for l in myset
        #             product *= pfmin[l,k]
        #         end # after this loop you got product==prod(pfmin[myset,:],dims=1)
        #         sumExp += log(1e-50 + product * prevminneg[j,k] + (1-prev[j,k]))
        #     end
        #     pfplus_[j] += ((-1)^length(myset)) * exp(sumExp)
        # end
    end
end; 
println("$t sec")

P_joint = pfplus[2:end,1] .* previn;
posterior = P_joint / pfplus[1];
println("Maximum difference pfplus: $(maximum(abs,pfplus_.-pfplus))")
println("Maximum difference posterior: $(maximum(abs,posterior_.-posterior))")



