include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

function get_patient_data()
    """
    This function stores the raw data from the Excel file from Alisa 
    """
    # Loading variables from Alisa's patient case (normal and high precision) and C++
    posterior_alisa, P_joint_alisa, dt_alisa = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    posterior_alisa_hp, P_joint_alisa_hp, dt_alisa_hp = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    posterior_cpp, P_joint_cpp, dt_cpp = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}(), Dict{String,Float64}()
    posterior_fortran, P_joint_fortran, dt_fortran = Dict{String,Vector{BigFloat}}(), Dict{String,Vector{BigFloat}}(), Dict{String,Float64}()
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
        # Load Fortran
        posterior_fortran["case $case_nr"] = readdlm("variables/fortran_patient404/output_case_$case_nr.csv",BigFloat)[1:end,1]
        P_joint_fortran["case $case_nr"] = readdlm("variables/fortran_patient404/output_case_$case_nr.csv",BigFloat)[1:end,2]
        dt_fortran["case $case_nr"] = readdlm("variables/fortran_patient404/output_case_$case_nr.csv")[1,3]
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

    methods_others = ["Matlab Sam","Matlab Alisa","Matlab Alisa hp","C++","Fortran"]
    posterior_others = [posterior_matlab,posterior_alisa,posterior_alisa_hp,posterior_cpp,posterior_fortran]
    P_joint_others = [P_joint_matlab,P_joint_alisa,P_joint_alisa_hp,P_joint_cpp,P_joint_fortran]
    dt_others = [dt_matlab,dt_alisa,dt_alisa_hp,dt_cpp,dt_fortran]
    return methods_others, posterior_others, P_joint_others, dt_others
end

function posteriors_func(case,methods,methods_others,posterior_others,P_joint_others,dt_others)
    """
    This function aims to process the raw data into one compact variable named 'posteriors' (containing more than Posteriors only)
    Then it runs the quickscore algorithm and stores these results too in this variable
    """
    # Retrieve the loaded data from the patient cases 
    patient_cases_raw, data_alisa, previn, pfmin, pfminneg = prepare_patient_data(case)
    # Save the data from other methods in the DataFrame
    posteriors = DataFrame(
        nr = collect(1:size(methods_others,1)), 
        Method = methods_others, 
        Posterior_min = Float64[minimum(posterior[case]) for posterior in posterior_others], 
        Posterior_max = Float64[maximum(posterior[case]) for posterior in posterior_others], 
        Posterior = Vector{BigFloat}[posterior[case] for posterior in posterior_others], 
        P_joint_min = Float64[minimum(P_joint[case]) for P_joint in P_joint_others], 
        P_joint_max = Float64[maximum(P_joint[case]) for P_joint in P_joint_others], 
        P_joint = Vector{Float64}[P_joint[case] for P_joint in P_joint_others], 
        time = Float64[dt[case] for dt in dt_others],
        previn = [Vector{Float64}([0.]) for _ in 1:size(methods_others,1)],
        pfmin = [Matrix{Float64}([[0.] [0.]]) for i=1:size(methods_others,1)],
        pfminneg = [Vector{Float64}([0.]) for _ in 1:size(methods_others,1)],
        pfplus = [Matrix{BigFloat}([[0.] [0.]]) for i=1:size(methods_others,1)]
    )
    # Run the quickscore algorithm for methods in Julia and store it in the DataFrame
    for (i,method) in enumerate(methods)
        pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg, method)
        push!(posteriors, (nr = i+size(method_others)[1],
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
methods_others, posterior_others, P_joint_others, dt_others = get_patient_data();
# Running the function for the different cases m=7,8,9
quickscore_methods = ["exp-sum-log","prod","prod QM Fl128","prod BF"]#,"exp-sum-log QM Fl128","prod BF Fl128"]#,"prod Fl32","exp-sum-log Fl32","prod BF Fl64","exp-sum-log BF"];
posteriors = Dict{String,DataFrame}(
    "case 1" => posteriors_func("case 1",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others),
    "case 2" => posteriors_func("case 2",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others),
    "case 3" => posteriors_func("case 3",quickscore_methods,method_others,posterior_others,P_joint_others,dt_others)
)

# Function that prints the data I want in nice tables
function print_posteriors(posteriors,case,m; lim1=1e-17, lim2=1e-6, save=false)    
    """
    This function prints the results of a certain patient case (m=7,8,9) using the variable 'posteriors' that contains all the data in a compact way.
    Via 'lim1' / 'lim2' you can set the coloring limits to highlight which quickscore methods give sensible results. 
    """
    posterior_diff = Float64.(maximum(abs,reshape(reduce(hcat,posteriors[case].Posterior),(size(posteriors[case].Posterior[1])[1],size(posteriors[case].Posterior)[1],1)) .- reshape(reduce(hcat,posteriors[case].Posterior),(size(posteriors[case].Posterior[1])[1],1,size(posteriors[case].Posterior)[1])),dims=1)[1,:,:])
    P_joint_diff = Float64.(maximum(abs,reshape(reduce(hcat,posteriors[case].P_joint),(size(posteriors[case].P_joint[1])[1],size(posteriors[case].P_joint)[1],1)) .- reshape(reduce(hcat,posteriors[case].P_joint),(size(posteriors[case].P_joint[1])[1],1,size(posteriors[case].P_joint)[1])),dims=1)[1,:,:])
    hl1 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] >= 0 && d[i, j] <= 1) : false, crayon"green")
    hl2 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (d[i, j] <= 0 || d[i, j] >= 1) : false, crayon"red")
    hl3 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim1 && d[i,j]!=0) : false, crayon"green")
    hl4 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) > lim1 && d[i,j]!=0) : false, crayon"red")
    hl5 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) < lim2 && d[i,j]!=0) : false, crayon"green")
    hl6 = Highlighter((d, i, j) -> typeof(d[i,j])==Float64 ? (abs(d[i, j]) > lim2 && d[i,j]!=0) : false, crayon"red")
    
    println("############################################## $case (m=$m) ##############################################")
    println("######################################################### P_joint & Posterior range ########################################################")
    pretty_table(posteriors[case][:,["nr","Method","P_joint_min","P_joint_max","Posterior_min","Posterior_max"]],alignment=:l,highlighters=(hl1,hl2))
    println("########################################################### P_joint (max-abs diff) ###########################################################")
    pretty_table(hcat(round.(P_joint_diff,sigdigits=2),string.(round.(posteriors[case].time,sigdigits=3))),header=vcat(posteriors[case].Method,"time"),row_names=posteriors[case].Method,alignment=:l,highlighters=(hl3,hl4))
    println("########################################################### posterior (max-abs diff) ###########################################################")
    pretty_table(hcat(round.(posterior_diff,sigdigits=5),string.(round.(posteriors[case].time,sigdigits=3))),header=vcat(posteriors[case].Method,"time"),row_names=posteriors[case].Method,alignment=:l,highlighters=(hl5,hl6))

    range = posteriors[case][:,["P_joint_min","P_joint_max","Posterior_min","Posterior_max","time"]];
    df_P_joint = DataFrame(hcat(posteriors[case].Method,P_joint_diff,posteriors[case].time),vcat("P_joint",posteriors[case].Method,"time"));
    df_posterior = DataFrame(hcat(posteriors[case].Method,posterior_diff,posteriors[case].time),vcat("posterior",posteriors[case].Method,"time"));

    if save
        CSV.write("variables/patient_404/methods_diff_case_$(parse(Int,case[end]))_range.csv",range)
        CSV.write("variables/patient_404/methods_diff_case_$(parse(Int,case[end]))_P_joint.csv",df_P_joint)
        CSV.write("variables/patient_404/methods_diff_case_$(parse(Int,case[end]))_posterior.csv",df_posterior)
        # XLSX.writetable("variables/patient_404/methods_diff_case_$(parse(Int,case[end])).xlsx","P_joint"=>df_P_joint,"posterior"=>df_posterior,"range"=>range; overwrite=true)
    end
    return range, df_P_joint, df_posterior
end 

# Calling the printing function for the different cases 
range_7, P_joint_diff_7, posterior_diff_7 = print_posteriors(posteriors,"case 1",7,save=false);
range_8, P_joint_diff_8, posterior_diff_8 = print_posteriors(posteriors,"case 2",8,save=false);
range_9, P_joint_diff_9, posterior_diff_9 = print_posteriors(posteriors,"case 3",9,save=false);


