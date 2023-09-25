using PrettyTables
using CSV
using DataFrames
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")
# previn,pfmin,pfminneg,actualdiseases,dt = quickscore_preparation(100,40,3,0.5)

# get stored variables from csv files 
file_vector = CSV.File("data_vectors.csv") 
file_matrix = CSV.File("data_matrix.csv") 
previn, pfminneg, actualdiseases = file_vector.previn, file_vector.pfminneg, file_vector.actualdiseases 
pfmin = Matrix(DataFrame(file_matrix)) 

pfplus_binary, posterior_binary, dt_binary = quickscore(previn, pfmin, pfminneg,true) 
pfplus_flip, posterior_flip, dt_flip = quickscore(previn, pfmin, pfminneg,false) 

# pretty_table(pfplus'/1e-8) 
pretty_table(hcat(previn,posterior_binary,actualdiseases),header=["prevalences", "posteriors", "actual diseases"])
pretty_table(hcat(previn,posterior_flip,actualdiseases),header=["prevalences", "posteriors", "actual diseases"])
print("Positive findings: $(size(pfmin,1)), time binary: $(dt_binary), time flip: $(dt_flip)")

