include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# LOAD SOME RANDOM VARIABLES previn,pfminneg,pfmin FOR A CERTAIN m=17 and run the quickscore
file = jldopen("variables/m_17.jld", "r")
    previn = file["previn"]
    pfminneg = file["pfminneg"]
    pfmin = file["pfmin"]
close(file)

m,n = size(pfmin)
prev = repeat(previn',inner=(n+1,1)); for i in 1:n prev[i+1,i]=1 end 
prevend = vcat(1,previn) # prevend are the prevalences 1, P(d_1=1),... P(d_n=1) and is needed to multiply at the end
prevminneg = prev.*pfminneg'

pfplus, P_joint, posterior,dt = quickscore(previn, pfmin, pfminneg,"prod");  





