include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")
include("useful_functions.jl")

previn, pfmin, pfminneg, actualdiseases, sens, sensneg = quickscore_preparation() 
sensneg 
prod(1 .- sensneg, dims=1)[1,:] 

matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
# patient = Dict()
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient_cases.mat"),"r")
    patient_matlab = read(matfile, "patient_to_julia")
    binary_cases = read(matfile, "binary_cases")[1,:]
    # diagtest_matlab = read(matfile, "diagtest_to_julia")
close(matfile)

index_patient = findall(sizeof.(patient_matlab) .== 0)
patient_matlab[index_patient] .= nothing
patient_array = dropdims(patient_matlab,dims=2)
patient_dict = Dict("Patient" => 1:size(patient_matlab,3), "age" => patient_array[1, :], "gender" => patient_array[2, :], "test" => patient_array[3, :], "tvalue" => patient_array[4, :], "diag" => patient_array[5, :], "dvalue" => patient_array[6, :], "cdiag" => patient_array[7, :], "ext" => patient_array[8, :], "evalue" => patient_array[9, :])
patient_df = DataFrame(patient_dict)
# patient_df[[306,307],[:Patient,:age,:cdiag]]
# patient_df[binary_cases,:]
patient_df[306,2]
patient_df[305:307,:]
patient = Patient[]

for i=1:size(patient_df,1)
    # print(i)
    age = patient_df[i,"age"] 
    gender = patient_df[i,"gender"]
    test = patient_df[i,"test"]
    tvalue = patient_df[i,"tvalue"]
    diag = patient_df[i,"diag"]
    dvalue = patient_df[i,"dvalue"]
    cdiag = patient_df[i,"cdiag"]
    ext = patient_df[i,"ext"]
    evalue = patient_df[i,"evalue"]
    println(typeof(age)," ",typeof(gender)," ",typeof(test)," ",typeof(tvalue)," ",typeof(diag)," ",typeof(dvalue)," ",typeof(cdiag)," ",typeof(evalue))
    Patient(age,gender,test,tvalue,diag,dvalue,cdiag,evalue)
    # push!(patient,Patient(age,gender,test,tvalue,diag,dvalue,cdiag,evalue))    
    # Patient(age,gender,test,tvalue,diag,dvalue,cdiag,evalue)
end


Patient(nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing,nothing)
# patient_df
patient_df[306,"evalue"]
# age = 70
# gender = 'M'
# test = [3284, 3812, 796, 3235, 1420]
# tvalue = ["true", "true", "true", "true", "true"]
# diag = Int[]  # An empty array for diag
# dvalue = String[]  # An empty array for dvalue
# cdiag = [2577, 0, 0, 0, 0, 0]
# ext = [2487, 1268]
# evalue = ["true", "true"]

patient_295 = 

patient_matlab


index_diagtest = findall(sizeof.(diagtest_matlab) .== 0)
diagtest_matlab[index_diagtest] .= nothing
# diagtest_dict = Dict("type" => spzeros(size(diagtest_matlab[1,:,:])), "mult" => spzeros(size(diagtest_matlab[2,:,:])), "sens" => spzeros(size(diagtest_matlab[3,:,:])))
idx = findall(diagtest_matlab[1,:,:] .!= nothing)
diagtest = spzeros(size(diagtest_matlab[1,:,:]))
for i in idx
    diagtest_dict[] = idx[i][1],idx[i][2]

end
# for (i,j)
# i
diagtest_dict["type"]
x[i] .= diagtest_matlab[1,i]
diagtest_matlab[1,i]
x[i]
i

x[i]
index2[1][1]

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

x[i]
diagtest_matlab[1,i]
x[i] .= diagtest_matlab[1,i]
diagtest_matlab


