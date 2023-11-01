include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# Loading data from a slightly adjusted version of Alisa's file 
matrix_alisa = Dict{String,Matrix{Any}}()
file = XLSX.openxlsx("/Users/sam/Downloads/patient cases.xlsx")
    matrix_alisa["case 1"] = file["Sheet1"][1:end, 1:end]
    matrix_alisa["case 2"] = file["Sheet2"][1:end, 1:end]
    matrix_alisa["case 3"] = file["Sheet3"][1:end, 1:end]
close(file)
# Loading saved variables from MATLAB
posterior_matlab, P_joint_matlab, dt_matlab = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/patient404.mat"),"r")
    for case_nr=1:3
        posterior_matlab["case $case_nr"] = read(matfile, "pdiag_case_$case_nr")[:,1]
        P_joint_matlab["case $case_nr"] = read(matfile, "P_joint_case_$case_nr")[:,1]
        dt_matlab["case $case_nr"] = read(matfile,"dt_case_$case_nr")
    end
close(matfile)



# Function to save the data in a dictionaries containing DataFrames 
function posteriors_func(case,methods,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)
    data_alisa = DataFrame(Float64.(matrix_alisa[case][8:end,:]),matrix_alisa[case][6,:])
    # diagnosis_nr = Int64.(data_alisa[:,"diagnosis ID number"])
    P_joint_hp_alisa = data_alisa[:,"Pjoint (high precision)"]
    P_joint_alisa = data_alisa[:,"Pjoint(standart double precision)"]
    posterior_hp_alisa = data_alisa[:,"Pdiag (high precision)"]
    posterior_alisa = data_alisa[:,"Pdiag(standart double precision)"]

    idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
    all_sens = Matrix(data_alisa[:,idx_sens])
    pos_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'true'"
    neg_idx = matrix_alisa[case][4,1:length(idx_sens)] .== "'false'"
    sens = all_sens[:,pos_idx]
    sensneg = all_sens[:,neg_idx]
    previn = data_alisa[:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]
    
    posteriors = DataFrame(
        nr = Int[1,2,3], 
        Method = String["Matlab Sam","Matlab Alisa","Matlab hp"], 
        # Posterior_range = Vector{Float64}[[minimum(posterior_matlab[case]),maximum(posterior_matlab[case])],[minimum(posterior_alisa),maximum(posterior_alisa)],[minimum(posterior_hp_alisa),maximum(posterior_hp_alisa)]], 
        Posterior_min = Float64[minimum(posterior_matlab[case]),minimum(posterior_alisa),minimum(posterior_hp_alisa)], 
        Posterior_max = Float64[maximum(posterior_matlab[case]),maximum(posterior_alisa),maximum(posterior_hp_alisa)], 
        Posterior = Vector{BigFloat}[posterior_matlab[case],posterior_alisa,posterior_hp_alisa], 
        # P_joint_range = Vector{Float64}[[minimum(P_joint_matlab[case]),maximum(P_joint_matlab[case])],[minimum(P_joint_alisa),maximum(P_joint_alisa)],[minimum(P_joint_hp_alisa),maximum(P_joint_hp_alisa)]], 
        P_joint_min = Float64[minimum(P_joint_matlab[case]),minimum(P_joint_alisa),minimum(P_joint_hp_alisa)], 
        P_joint_max = Float64[maximum(P_joint_matlab[case]),maximum(P_joint_alisa),maximum(P_joint_hp_alisa)], 
        P_joint = Vector{BigFloat}[P_joint_matlab[case],P_joint_alisa,P_joint_hp_alisa],
        time = Float64[dt_matlab[case],matrix_alisa[case][7,4],matrix_alisa[case][7,3]],
        previn = Vector{Float64}[[0.],[0.],[0.]],
        pfmin = Matrix{Float64}[[[0.] [0.]],[[0.] [0.]],[[0.] [0.]]],
        pfminneg = Vector{Float64}[[0.],[0.],[0.]],
        pfplus = Matrix{BigFloat}[[[0.] [0.]],[[0.] [0.]],[[0.] [0.]]]
    )

    for (i,method) in enumerate(methods)
        println(method)
        pfplus, posterior, dt = quickscore(previn, pfmin, pfminneg, true, method)
        push!(posteriors, (nr = i+3,
            Method = method, 
            # Posterior_range = [minimum(posterior),maximum(posterior)],
            Posterior_min = minimum(posterior),
            Posterior_max = maximum(posterior),
            Posterior = posterior, 
            # P_joint_range = [minimum(pfplus[2:end,1].*previn),maximum(pfplus[2:end,1].*previn)], 
            P_joint_min = minimum(pfplus[2:end,1].*previn),
            P_joint_max = maximum(pfplus[2:end,1].*previn),
            P_joint = pfplus[2:end,1].*previn,
            time = dt,
            previn = previn,
            pfmin = pfmin,
            pfminneg = pfminneg,
            pfplus = pfplus))
    end

    return posteriors
end

# Running the function for the different cases m=7,8,9
methods = ["exp-sum-log","exp-sum-log BF(all)","prod","prod BF(all)","prod BF(pfmin)"]#,"exp-sum-log BF(pfmin)"]
posteriors = Dict{String,DataFrame}(
    "case 1" => posteriors_func("case 1",methods,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab),
    "case 2" => posteriors_func("case 2",methods,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab),
    "case 3" => posteriors_func("case 3",methods,matrix_alisa,P_joint_matlab,posterior_matlab,dt_matlab)
)




# Function that prints the data I want in nice tables
function print_posteriors(posteriors,case,m,lim1,lim2)    
    diff_posterior, diff_P_joint = Vector{Float64}[], Vector{Float64}[]
    for i=eachindex(posteriors[case].Method)
        push!(diff_posterior,Float64.(maximum(abs,reduce(hcat,posteriors[case].Posterior) .- posteriors[case].Posterior[i],dims=1)[1,:]))
        push!(diff_P_joint,Float64.(maximum(abs,reduce(hcat,posteriors[case].P_joint) .- posteriors[case].P_joint[i],dims=1)[1,:]))
    end
    hl1 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] >= 0 && d[i, j] <= 1) : false, crayon"green")
    hl2 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] <= 0 || d[i, j] >= 1) : false, crayon"red")
    hl3 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim1 && d[i,j]!=0) : false, crayon"green")
    hl4 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim2 && d[i,j]!=0) : false, crayon"green")

    println("############################################## $case (m=$m) ##############################################")
    println("############################################## range of P_joint and Posterior between methods ##############################################")
    pretty_table(posteriors[case][:,["nr","Method","P_joint_min","P_joint_max","Posterior_min","Posterior_max"]],alignment=:l,highlighters=(hl1,hl2))
    println("############################################## max-abs diffference of P_joint between methods ##############################################")
    pretty_table(hcat(reduce(hcat,diff_P_joint),posteriors[case].time),row_names=posteriors[case].Method,header=vcat(posteriors[case].Method,"time"),alignment=:l,highlighters=hl3) # ,formatters=ft_printf("%1.0e")
    println("########################################### max-abs diffference of Posterior between methods ###########################################")
    pretty_table(hcat(reduce(hcat,diff_posterior),posteriors[case].time),row_names=posteriors[case].Method,header=vcat(posteriors[case].Method,"time"),alignment=:l,highlighters=hl4)
end 

# Calling the printing function for the different cases 
print_posteriors(posteriors,"case 1",7,1e-17,1e-6)
print_posteriors(posteriors,"case 2",8,1e-17,1e-6)
print_posteriors(posteriors,"case 3",9,1e-17,1e-6)


# This shows that it does not matter what precision in BigFloat you choose, it is always slow 
pfmin = posteriors["case 1"].pfmin[4] # pick pfmin of 'case 1'
println("normal BigFloat")
@benchmark prod(BigFloat.(pfmin),dims=1);
println("normal BigFloat with specification of precision = 256 bits")
@btime prod(BigFloat.(pfmin,precision=256),dims=1);
println("precision of 100 bits")
@btime prod(BigFloat.(pfmin,precision=100),dims=1);
println("same precision = 53 bits as Float64")
@btime prod(BigFloat.(pfmin,precision=53),dims=1);
println("normal Float64 (no BigFloat at all)")
@btime prod(pfmin,dims=1);

# # This shows that posterior=pfplus[2:end,1].*previn./pfplus[1] is still bad if you take only the denominator pfplus[1] using BigFloat(pfmin)
# p1 = posteriors["case 1"].pfplus[4][2:end,1] .* posteriors["case 1"].previn[4] ./ posteriors["case 1"].pfplus[8][1];
# p2 = posteriors["case 1"].pfplus[6][2:end,1] .* posteriors["case 1"].previn[6] ./ posteriors["case 1"].pfplus[8][1];
# p3 = posteriors["case 1"].pfplus[8][2:end,1] .* posteriors["case 1"].previn[8] ./ posteriors["case 1"].pfplus[8][1];
# p4 = posteriors["case 1"].Posterior[8];
# println("Difference\t 'prod BF(pfmin)'\t and\t 'exp-sum-log + BF(pfmin)@only pfplus[1]':\t $(Float64(maximum(abs,p4 - p1)))")
# println("Difference\t 'prod BF(pfmin)'\t and\t 'prod + BF(pfmin)_for_pfplus[1]':\t\t $(Float64(maximum(abs,p3 - p2)))")
# println("Difference\t 'prod BF(pfmin)'\t and\t 'prod BF(pfmin)' (calculated differently):\t $(maximum(abs,p4 - p3))")



