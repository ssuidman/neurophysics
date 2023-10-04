include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# # THIS IS A WAY TO RUN THE CODE FOR AN EXAMPLE 
# m=16
# file = jldopen("variables/m_$m.jld", "r")
# previn = file["previn"]
# pfminneg = file["pfminneg"]
# pfmin = file["pfmin"]
# close(file)
# pfplus,posterior,dt = quickscore(previn, pfmin, pfminneg, false)


# THESE ARE SOME PROFILING TOOLS
# @pprof pfplus, posterior, dt = quickscore(previn,pfmin,pfminneg,true)
# @code_warntype quickscore(previn, pfmin, pfminneg, false)
# @benchmark log.(rand(101,100))


# # THIS LINE PRINTS OUT WHAT IT IS PRECOMPILING WHEN RUNNING IN THE TERMINAL 
# julia --trace-compile=stderr main.jl 


# THIS IS A COMPARISON BETWEEN MOST IMPORTANT LINES 
function f()
    x, myset, pfmin, prevminneg, prev = rand(101, 100), [1,4,5,6], rand(18,100), rand(101,100), rand(101,100)

    t = @elapsed begin
        for i = 1:2^24
            # y = log.(x)
            y = prod(x,dims=1)

            # y = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
            # y = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
            # y = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))
            # y = (1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))
            # y = prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
        end
    end
    println("Running time was $t seconds")
end

f()

# For taking log of matrix(101,100)         [test 2^20       Julia.ipynb: 42.53 sec      Matlab: 24.87 sec     Julia.jl: 41.58]
# For taking prod of matrix(101,100)        [test 2^24       Julia.ipynb: 24.23 sec      Matlab: 74.43 sec     Julia.jl: 24.71] 

# Quickscore total log                      [test 2^20       Julia.ipynb: 62.19 sec      Matlab: 39.98 sec     Julia.jl: 61.78]
# Quickscore total prod                     [test 2^22       Julia.ipynb: 49.78 sec      Matlab: 47.48 sec     Julia.jl: 44.88]
# Quickscore including taking the log       [test 2^20       Julia.ipynb: 58.45 sec      Matlab: 33.69 sec     Julia.jl: 62.51]
# Quickscore excluding taking the log       [test 2^22       Julia.ipynb: 37.33 sec      Matlab: 35.13 sec     Julia.jl: 39.13]
# Quickscore taking the prod                [test 2^22       Julia.ipynb: 47.17 sec      Matlab: 47.11 sec     Julia.jl: 48.18]


