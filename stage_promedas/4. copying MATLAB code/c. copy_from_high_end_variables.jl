include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../0. useful/useful_functions.jl")
include("b. copying_functions.jl")

patient, test, diagn, diagtest, prev, belowSens, aboveSens, restSens, malemult, femalemult, agemult, parents, link_ext = import_matlab_variables(); 
# diagtest[3421,919] #[3216,3997] [1634,3325] [2351,3470]
# agemult[77:81,:]
# i = 295
i = 404 
patient404 = DataFrame(patient[404,:])
patient404[1,"test"] = [2867,2502,1731,3568,437,3426,3978,1122,1019,55,3118,987]
patient404[1,"tvalue"] = ["true","true","true","false","true","true","false","true","true","false","false","false"]
diagset404 = sort(unique(reduce(vcat,parents[patient404[1,"test"][patient404[1,"tvalue"].=="true"]])))
belowSens
aboveSens
restSens
# diagset = sort(unique(reduce(vcat,parents[patient[i,"test"]])))
size(diagset404)
diagset2_matlab
diagset404

setdiff(union(diagset2_matlab,diagset404),intersect(diagset2_matlab,diagset404))
1228 in parents[patient404[1,"test"]][6]
patient404[1,"test"][6]
test[3426,:]

bit = zeros(4011) .== 1
bit[diagset404] .= 1
sum(test.type .== 0 .&& bit)
size(x)
size(diagset404)


patientprev = patient[i,"gender"] !== nothing ? (patient[i,"gender"] == "M" ? prev.*malemult : prev.*femalemult) : prev 
for (j,d) in enumerate(eachrow(agemult))
    if typeof(d[1]) == Float64
        mult = (d[1] .<= patient[i,"age"]) .&& (patient[i,"age"] .< d[2]) ? d[3] : 1
    elseif typeof(d[1]) == Vector{Float64}
        mult_idx = (d[1] .<= patient[i,"age"]) .&& (patient[i,"age"] .< d[2]) 
        mult = sizeof(d[3][mult_idx]) != 0 ? d[3][mult_idx][1] : 1
    else
        mult = 1
    end
    patientprev[j] *= mult
end # d = 3612, 3613, 78, 79, 3566

if patient[i,"ext"] != nothing
    test_ext = patient[i,"ext"][patient[i,"evalue"] .== "true"]
    for diag_ext in test_ext
        idx = filter(ext_idx -> ext_idx[2] in diag_ext, findall(!iszero,link_ext))
        d = [k[1] for k in idx]
        mult = DataFrame(diagtest[idx]).mult
        patientprev[d] .*= mult
    end
end 


patientprev[diagset404]
patientprev[diagset2_matlab]

# diagset404[patientprev[diagset404] .!=0]

# CartesianIndex(36, 2487)
# CartesianIndex(149, 2487)
# CartesianIndex(174, 2487)
# CartesianIndex(262, 2487)
# CartesianIndex(3637, 1268)
# CartesianIndex(3644, 1268)
# CartesianIndex(3658, 1268)
# CartesianIndex(3683, 1268)











