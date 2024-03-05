include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")


# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick Fl128");  

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

