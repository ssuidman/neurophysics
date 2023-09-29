using PrettyTables
using CSV
using DataFrames
using JLD2
using MAT
using Plots
using Profile
using PProf
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

file = jldopen("variables/m_15.jld", "r")
previn = file["previn"]
pfminneg = file["pfminneg"]
pfmin = file["pfmin"]
close(file)

# This gives a link to the scheme that shows which parts cost most time
# Sometimes you need to (restart and) run everything from scratch to show the link 
pfplus,posterior,dt = 0,0,0
@pprof pfplus, posterior, dt = quickscore(previn,pfmin,pfminneg,true)
