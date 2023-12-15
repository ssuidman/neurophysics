include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")

# posterior_diff = []
# for i=1:20
#     previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(i,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
#     pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick BF thread"); 
#     pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread"); 
#     push!(posterior_diff,maximum(abs,posterior .- posterior_))
# end

# # Write variables to jld-file 
# jldopen("variables/error_BF_Fl128.jld","w") do file
    # file["posterior_diff"] = posterior_diff
    # file["pfplus_diff"] = pfplus_diff
# end

# Retrieve variables from jld-file 
jldopen("variables/error_BF_Fl128.jld", "r") do file
    global posterior_diff = file["posterior_diff"]
    # global pfplus_diff = file["pfplus_diff"]
end

plot(1:20,posterior_diff,yscale=:log10)
title!("max-abs-diff across diagnoses for BF & Fl128")
xlabel!("m (positive tests)")
ylabel!("diff")

