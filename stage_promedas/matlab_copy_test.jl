include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

matlab_path = "/Users/sam/Documents/MATLAB/code/"

# read pumcu_diagsIdName.txt as Vector{String} sorted along first dimension
pumcu_diagsIdName = readdlm(matlab_path*"db79/pumcu_diagsIdName.txt",'@') 
diagn = map(string,sortslices(pumcu_diagsIdName,dims=1)[:,3])
# read pumcu_testsIdName.txt
pumcu_testsIdName = readdlm(matlab_path*"db79/pumcu_testsIdName.txt",'@')
test = Dict{Any, Vector{}}()
test["gender"] = map(string,sortslices(pumcu_testsIdName,dims=1)[:,2])
test["name"] = map(string,sortslices(pumcu_testsIdName,dims=1)[:,3])
# read pumcu_externalsIdName.txt 
pumcu_externalsIdName = readdlm(matlab_path*"db79/pumcu_externalsIdName.txt",'@')
test["type"] = zeros(Int64,size(pumcu_testsIdName,1))
test["type"][pumcu_externalsIdName[:,1]] .= 6
# read pumcu_diagprevs.txt
pumcu_diagprevs = readdlm(matlab_path*"db79/pumcu_diagprevs.txt",'@')
prev = sortslices(pumcu_diagprevs,dims=1)[:,2]
# read pumcu_agemult.txt
pumcu_agemult = readdlm(matlab_path*"db79/pumcu_agemult.txt",' ')
agemult_matrix = sortslices(pumcu_agemult,dims=1)

x = sort(collect(countmap(agemult_matrix[:,1])))
for (i,d) in enumerate(x)
    if x[i][2] > 1
        double_index = findall(agemult_matrix[:,1].==x[i][1])
        println(typeof(agemult_matrix[double_index,:])," ",double_index,"--> ",agemult_matrix[double_index,:])
    end
end

double_index = findall(agemult_matrix[:,1].==78.0)
for i in double_index
    println(agemult_matrix[i,:])
end
