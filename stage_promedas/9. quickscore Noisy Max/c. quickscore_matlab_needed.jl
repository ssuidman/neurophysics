include("../0. useful/quickscore_algorithm.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"MATLAB"); 

matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"; 
matfile = matopen(matlab_dir*"variables/noisy_max/patient_live_input.mat","r"); 
    patientprev = read(matfile, "patientprev")[:,1]; 
    sens_normal = read(matfile, "sens_normal")[:,1]; 
    sens_medium = read(matfile, "sens_medium"); 
    sens_extreme = read(matfile, "sens_extreme"); 
close(matfile)

previn = copy(patientprev); 
pfmin = vcat(prod(sens_medium,dims=3)[:,:,1]',sens_extreme'); 
pfminneg = copy(sens_normal); 
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"Fl128 trick thread"); 

matfile = matopen(matlab_dir*"variables/noisy_max/patient_live_output.mat","w")
    write(matfile,"pdiag",Float64.(posterior))
close(matfile)


