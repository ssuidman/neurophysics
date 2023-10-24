include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")
include("useful_functions.jl")

# previn, pfmin, pfminneg, actualdiseases, sens, sensneg = quickscore_preparation() 
# sensneg 
# prod(1 .- sensneg, dims=1)[1,:] 


function patient_func(patient_matlab::Array{Any,3})
    patient_dict = Dict{String, Vector{Any}}()
    for i in axes(patient_matlab,3) 
        age = sizeof(patient_matlab[1,1,i]) != 0 ? Int(patient_matlab[1,1,i]) : nothing
        gender = sizeof(patient_matlab[2,1,i]) != 0 ? String(patient_matlab[2,1,i]) : nothing
        test = sizeof(patient_matlab[3,1,i]) != 0 ? (typeof(patient_matlab[3,1,i]) == Matrix{Float64} ? Int.(patient_matlab[3,1,i][1,:]) : Int(patient_matlab[3,1,i])) : nothing
        if !any(x -> typeof(x) == Matrix{Union{}}, patient_matlab[4,1,i])
            tvalue = sizeof(patient_matlab[4,1,i]) != 0 ? String.(patient_matlab[4,1,i][1,:]) : nothing
        else 
            tvalue = String.(filter(item -> !(typeof(item) == Matrix{Union{}}), patient_matlab[4,1,i]))
        end
        diag = sizeof(patient_matlab[5,1,i]) != 0 ? (typeof(patient_matlab[5,1,i]) == Matrix{Float64} ? Int.(patient_matlab[5,1,i][1,:]) : Int(patient_matlab[5,1,i])) : nothing
        dvalue = sizeof(patient_matlab[6,1,i]) != 0 ? String.(patient_matlab[6,1,i][1,:]) : nothing
        cdiag = sizeof(patient_matlab[7,1,i]) != 0 ? (typeof(patient_matlab[7,1,i]) == Matrix{Float64} ? Int.(patient_matlab[7,1,i][1,:]) : Int(patient_matlab[7,1,i])) : nothing
        ext = sizeof(patient_matlab[8,1,i]) != 0 ? (typeof(patient_matlab[8,1,i]) == Matrix{Float64} ? Int.(patient_matlab[8,1,i][1,:]) : Int(patient_matlab[8,1,i])) : nothing
        evalue = sizeof(patient_matlab[9,1,i]) != 0 ? String.(patient_matlab[9,1,i][1,:]) : nothing
        patient_dict["age"] = push!(get(patient_dict, "age", []), age)
        patient_dict["gender"] = push!(get(patient_dict, "gender", []), gender)
        patient_dict["test"] = push!(get(patient_dict, "test", []), test)
        patient_dict["tvalue"] = push!(get(patient_dict, "tvalue", []), tvalue)
        patient_dict["diag"] = push!(get(patient_dict, "diag", []), diag)
        patient_dict["dvalue"] = push!(get(patient_dict, "dvalue", []), dvalue)
        patient_dict["cdiag"] = push!(get(patient_dict, "cdiag", []), cdiag)
        patient_dict["ext"] = push!(get(patient_dict, "ext", []), ext)
        patient_dict["evalue"] = push!(get(patient_dict, "evalue", []), evalue)
    end
    patient = DataFrame(patient_dict)
    return patient
end



function diagtest_func(diagtest_matlab::Array{Any,3})
    idx = findall(sizeof.(diagtest_matlab[1,:,:]) .!= 0)
    diagtest = Array{Any}(nothing, size(diagtest_matlab,2), size(diagtest_matlab,3))
    for i in idx
        type = typeof(diagtest_matlab[1,i]) == Float64 ? [Int(diagtest_matlab[1,i])] : Int.(diagtest_matlab[1,i][1,:]) # gefikst
        mult = sizeof(diagtest_matlab[2,i]) != 0 ? diagtest_matlab[2,i][1] : nothing
        sens = sizeof(diagtest_matlab[3,i]) != 0 ? Vector(diagtest_matlab[3,i][1,:]) : nothing
        diagtest[i] = Dict("type"=>type, "mult"=>mult, "sens"=>sens)
    end
    return diagtest
end



function test_func(test_matlab::Array{Any,3})
    test_dict = Dict{String, Vector{Any}}()
    for i in axes(test_matlab,3) 
        gender = test_matlab[1,1,i]
        name = test_matlab[2,1,i]
        type = typeof(test_matlab[3,1,i]) == Float64 ? test_matlab[3,1,i] : test_matlab[3,1,i][1,:]
        normal = sizeof(test_matlab[4,1,i]) != 0 ? test_matlab[4,1,i][1,:] : nothing
        test_dict["gender"] = push!(get(test_dict, "gender", []), gender)
        test_dict["name"] = push!(get(test_dict, "name", []), name)
        test_dict["type"] = push!(get(test_dict, "type", []), type)
        test_dict["normal"] = push!(get(test_dict, "normal", []), normal)
    end
    test = DataFrame(test_dict)
    return test
end



function agemult_func(agemult_matlab::Array{Any,3})
    agemult_array = copy(agemult_matlab[:,1,:]')
    idx1 = sizeof.(agemult_matlab) .== 0
    agemult_array[idx1[:,1,:]'] .= nothing
    idx2 = (typeof.(agemult_array) .!= Float64) .&& (typeof.(agemult_array) .!= Nothing)
    agemult_array[idx2] .= Vector{Vector{Float64}}(vec.(agemult_array[idx2]))
    agemult = DataFrame(agemult_array,["low","high","mult"])
    return agemult
end



function parents_func(parents_matlab::Vector{Any})
    parents = Vector{Any}(undef, 4011)
    for i in axes(parents_matlab,1)
        parents[i] = sizeof(parents_matlab[i]) != 0 ? (typeof(parents_matlab[i]) == Matrix{Float64} ? Int.(parents_matlab[i][:,1]) : Int(parents_matlab[i])) : nothing
    end
    return parents
end



function import_matlab_variables()
    matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
    matfile = matopen(joinpath(matlab_dir,"variables_julia/patient_cases.mat"),"r")
        patient_matlab = read(matfile, "patient_to_julia")
        test_matlab = read(matfile, "test_to_julia")
        diagn = read(matfile,"diagn")[1,:]
        diagtest_matlab = read(matfile, "diagtest_to_julia")
        prev = read(matfile, "prev")[1,:]
        belowSens = read(matfile, "belowSens")
        aboveSens = read(matfile, "aboveSens")
        restSens = read(matfile, "restSens")
        malemult = read(matfile, "malemult")[1,:]
        femalemult = read(matfile, "femalemult")[1,:]
        agemult_matlab = read(matfile, "agemult_to_julia")
        parents_matlab = read(matfile, "parents")[1,:]
        link_ext = read(matfile, "link_ext")
    close(matfile)

    patient = patient_func(patient_matlab)
    diagtest = diagtest_func(diagtest_matlab)
    test = test_func(test_matlab)
    agemult = agemult_func(agemult_matlab)
    parents = parents_func(parents_matlab)
    
    return patient, test, diagn, diagtest, prev, belowSens, aboveSens, restSens, malemult, femalemult, agemult, parents, link_ext
end

