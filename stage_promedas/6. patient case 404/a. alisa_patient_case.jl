include("../packages.jl")
include("../quickscore_algorithm.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

function get_patient_data()
    """
    This function stores the raw data from the Excel file from Alisa 
    """
    # Loading variables from Alisa's patient case (normal and high precision) and C++
    posterior_alisa, P_joint_alisa, dt_alisa = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    posterior_alisa_hp, P_joint_alisa_hp, dt_alisa_hp = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    posterior_cpp, P_joint_cpp, dt_cpp = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    for case_nr=1:3
        # Get the raw data from the function in the C++ folder
        patient_cases_raw, data_alisa, previn, pfmin, pfminneg = prepare_patient_data("case $case_nr")
        # Load Alisa normal 
        posterior_alisa["case $case_nr"] = data_alisa[:,"Pdiag(standart double precision)"]
        P_joint_alisa["case $case_nr"] = data_alisa[:,"Pjoint(standart double precision)"]
        dt_alisa["case $case_nr"] = patient_cases_raw["case $case_nr"][7,4]
        # Load Alisa hp 
        posterior_alisa_hp["case $case_nr"] = data_alisa[:,"Pdiag (high precision)"]
        P_joint_alisa_hp["case $case_nr"] = data_alisa[:,"Pjoint (high precision)"]
        dt_alisa_hp["case $case_nr"] = patient_cases_raw["case $case_nr"][7,3]
        # Load C++ 
        posterior_cpp["case $case_nr"] = Float64.(readdlm("variables/cpp/cpp_output_case_$case_nr.csv",',')[2:end,1])
        P_joint_cpp["case $case_nr"] = Float64.(readdlm("variables/cpp/cpp_output_case_$case_nr.csv",',')[2:end,2])
        dt_cpp["case $case_nr"] = readdlm("variables/cpp/cpp_output_case_$case_nr.csv",',')[2,3]
    end
    # Load my own MATLAB run 
    posterior_matlab, P_joint_matlab, dt_matlab = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
    matfile = matopen(joinpath(matlab_dir,"variables/patient_404/patient404.mat"),"r")
        for case_nr=1:3
            posterior_matlab["case $case_nr"] = read(matfile, "pdiag_case_$case_nr")[:,1]
            P_joint_matlab["case $case_nr"] = read(matfile, "P_joint_case_$case_nr")[:,1]
            dt_matlab["case $case_nr"] = read(matfile,"dt_case_$case_nr")
        end
    close(matfile)

    method_others = ["Matlab Sam","Matlab Alisa","Matlab Alisa hp","C++"]
    posterior_others = [posterior_matlab,posterior_alisa,posterior_alisa_hp,posterior_cpp]
    P_joint_others = [P_joint_matlab,P_joint_alisa,P_joint_alisa_hp,P_joint_cpp]
    dt_others = [dt_matlab,dt_alisa,dt_alisa_hp,dt_cpp]
    return method_others, posterior_others, P_joint_others, dt_others
end



function posteriors_func(case,methods,method_others,posterior_others,P_joint_others,dt_others)
    """
    This function aims to process the raw data into one compact variable named 'posteriors' (containing more than Posteriors only)
    Then it runs the quickscore algorithm and stores these results too in this variable
    """
    # Retrieve the loaded data from the patient cases 
    patient_cases_raw, data_alisa, previn, pfmin, pfminneg = prepare_patient_data(case)
    # Save the data from other methods in the DataFrame
    posteriors = DataFrame(
        nr = Int[1,2,3,4], 
        Method = method_others, 
        Posterior_min = Float64[minimum(posterior[case]) for posterior in posterior_others], 
        Posterior_max = Float64[maximum(posterior[case]) for posterior in posterior_others], 
        Posterior = Vector{BigFloat}[posterior[case] for posterior in posterior_others], 
        P_joint_min = Float64[minimum(P_joint[case]) for P_joint in P_joint_others], 
        P_joint_max = Float64[maximum(P_joint[case]) for P_joint in P_joint_others], 
        P_joint = Vector{Float64}[P_joint[case] for P_joint in P_joint_others], 
        time = Float64[dt[case] for dt in dt_others],
        previn = Vector{Float64}[[0.],[0.],[0.],[0.]],
        pfmin = Matrix{Float64}[[[0.] [0.]],[[0.] [0.]],[[0.] [0.]],[[0.] [0.]]],
        pfminneg = Vector{Float64}[[0.],[0.],[0.],[0.]],
        pfplus = Matrix{BigFloat}[[[0.] [0.]],[[0.] [0.]],[[0.] [0.]],[[0.] [0.]]]
    )
    # Run the quickscore algorithm for methods in Julia and store it in the DataFrame
    for (i,method) in enumerate(methods)
        pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg, method)
        push!(posteriors, (nr = i+3,
            Method = method, 
            Posterior_min = minimum(posterior),
            Posterior_max = maximum(posterior),
            Posterior = posterior, 
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

# Retrieve the data from the patient cases and others
method_others, posterior_others, P_joint_others, dt_others = get_patient_data()
# Running the function for the different cases m=7,8,9
quickscore_methods = ["exp-sum-log","prod","prod Fl32","exp-sum-log Fl32","prod BF","prod BF Fl128"]#,"prod BF Fl64","exp-sum-log BF"]
posteriors = Dict{String,DataFrame}(
    "case 1" => posteriors_func("case 1",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others),
    "case 2" => posteriors_func("case 2",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others),
    "case 3" => posteriors_func("case 3",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others)
)

# Function that prints the data I want in nice tables
function print_posteriors(posteriors,case,m,lim1,lim2)    
    """
    This function prints the results of a certain patient case (m=7,8,9) using the variable 'posteriors' that contains all the data in a compact way.
    Via 'lim1' / 'lim2' you can set the coloring limits to highlight which quickscore methods give sensible results. 
    """
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





