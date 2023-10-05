include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")



# # STORING THE VARIABLES IN JL-FILES AND MAT-FILES, DO NOT RUN THIS ONE OFTEN, BECAUSE PREVIN,PFMINNEG,PFMIN THEN CHANGE!!!
# for i in 1:22
#     previn, pfmin, pfminneg, actualdiseases = quickscore_preparation(100,30,3,0.5,i)
    
#     #write variables to jld-file 
#     jldopen("variables/m_$i.jld","w") do file
#         file["previn"] = previn
#         file["pfminneg"] = pfminneg
#         file["actualdiseases"] = actualdiseases
#         file["pfmin"] = pfmin
#     end
    
#     # Write variables to mat-file
#     matfile = matopen("/Users/sam/Documents/MATLAB/Wim/variables_julia/m_$i.mat","w")
#     write(matfile,"previn",previn)
#     write(matfile,"pfminneg",pfminneg)
#     write(matfile,"actualdiseases",actualdiseases)
#     write(matfile,"pfmin",pfmin)
#     close(matfile)
# end



# # RUNNING THE QUICKSCORE ALGORITHM FROM JL-FILES FOR ALL M=1,...,22
# for i in 1:22
#     jldopen("variables/m_$i.jld", "r") do file
#         global previn = file["previn"]
#         global pfminneg = file["pfminneg"] 
#         global actualdiseases = file["actualdiseases"] 
#         global pfmin = file["pfmin"] 
#     end
    
#     pfplus,posterior,dt = quickscore(previn, pfmin, pfminneg,false) 

#     jldopen("variables/m_$i.jld","a") do file # make sure that this part has not yet run, otherwise it gives an error that parts of the file already exist
#         file["pfplus"] = pfplus 
#         file["posterior"] = posterior
#         file["dt"] = dt
#     end
# end



# LOAD SOME RANDOM VARIABLES previn,pfminneg,pfmin FOR A CERTAIN m=17 and run the quickscore
# file = jldopen("variables/m_17.jld", "r")
# previn = file["previn"]
# pfminneg = file["pfminneg"]
# pfmin = file["pfmin"]
# close(file)
# pfplus,posterior,dt = quickscore(previn, pfmin, pfminneg,false) 


# RUN THE MATLAB QUICKSCORE ALGORITHM IN THE FILE "run_quickscore_1_22.m"



# # GETTING THE VARIABLES FROM THE JLD-FILES AND MAT-FILES TO BE ABLE TO WORK WITH THEM --> START HERE NEXT TIME 
# using JLD2
# using MAT
# dt_julia = Vector{Float64}()
# dt_matlab = Vector{Float64}()
# posterior_julia = Vector{Vector{Float64}}()
# posterior_matlab = Vector{Vector{Float64}}()

# for i in 1:22
#     jlfile = jldopen("variables/m_$i.jld", "r")
#     push!(dt_julia, jlfile["dt"])
#     posterior = jlfile["posterior"]
#     push!(posterior_julia,posterior)
#     close(jlfile)
    
#     matfile = matopen("/Users/sam/Documents/MATLAB/Wim/variables_julia/m_$i.mat", "r")
#     dt = read(matfile, "timep")
#     push!(dt_matlab,dt)
#     posterior = vec(read(matfile, "posterior"))
#     push!(posterior_matlab,posterior)
#     close(matfile)
# end
# posterior_julia = hcat(posterior_julia...)
# posterior_matlab = hcat(posterior_matlab...)



# # MAKE PLOTS 
# using Plots
# gr()
# plot!()

# p1 = plot(1:22,dt_julia,label="Julia",legend=:topleft)
# plot!(1:22,dt_matlab,label="MATLAB")
# xlabel!("postive tests")
# ylabel!("t(s)")

# p2 = plot(1:22,dt_julia,label="Julia",yscale= :log10,legend=:topleft)
# plot!(1:22,dt_matlab,label="MATLAB",yscale= :log10)
# xlabel!("postive tests")

# plot(p1,p2,layout=(1,2),suptitle="Running time quickscore")
# # savefig("figures/running_time_comparison_quickscore.png")







