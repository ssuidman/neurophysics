using CSV
using DataFrames

data_vector = DataFrame(previn=previn,pfminneg=pfminneg,actualdiseases=actualdiseases)
data_matrix = DataFrame(pfmin, :auto)
CSV.write("data_vectors.csv",data_vector)
CSV.write("data_matrix.csv",data_matrix)

file_vector = CSV.File("data_vectors.csv")
file_matrix = CSV.File("data_matrix.csv")
previn = file_vector.previn
pfminneg = file_vector.pfminneg
actualdiseases = file_vector.actualdiseases
pfmin = Matrix(DataFrame(file_matrix))
# previn = file.previn