include("packages.jl")
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"

# Write variables to jld-file 
jldopen("variables/data.jld","w") do file
    file["previn"] = previn
    file["pfminneg"] = pfminneg
    file["actualdiseases"] = actualdiseases
    file["pfmin"] = pfmin
end

# Retrieve variables from jld-file 
jldopen("variables/data.jld", "r") do file
    global previn = file["previn"]
    global pfminneg = file["pfminneg"] 
    global actualdiseases = file["actualdiseases"] 
    global pfmin = file["pfmin"] 
end

# Write variables to mat-file
matfile = matopen(joinpath(matlab_dir,"variables/test.mat"),"w")
    write(matfile,"previn",previn)
    write(matfile,"pfminneg",pfminneg)
    write(matfile,"actualdiseases",actualdiseases)
    write(matfile,"pfmin",pfmin)
close(matfile)

# Retrive variables from mat-file 
matfile = matopen(joinpath(matlab_dir,"variables/test.mat"),"r")
    previn = read(matfile, "previn")
    pfminneg = read(matfile, "pfminneg")
    actualdiseases = read(matfile, "actualdiseases")
    pfmin = read(matfile, "pfmin")
close(matfile)




