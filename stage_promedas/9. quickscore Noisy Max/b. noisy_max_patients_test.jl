include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../4. copying MATLAB code/b. copying_functions.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"MATLAB"); 

cases = ["128","404_case_1","404_case_2","404_case_3","6","42"]
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"; 
matfile = matopen(joinpath(matlab_dir,"variables/noisy_max/patient_$(cases[6]).mat"),"r"); 
    patientprev = read(matfile, "patientprev")[:,1]; 
    sens_normal = read(matfile, "sens_normal")[:,1]; 
    sens_medium = read(matfile, "sens_medium"); 
    sens_extreme = read(matfile, "sens_extreme"); 
    diagset = Int.(read(matfile, "diagset"))[:,1]; 
    pdiag = read(matfile, "pdiag")[:,1]; 
close(matfile)

pfmin = vcat(prod(sens_medium,dims=3)[:,:,1]',sens_extreme')
pfminneg = copy(sens_normal)
previn = copy(patientprev)
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"Fl128 trick thread");

file = XLSX.openxlsx("/Users/sam/Large_files/Data/stage_promedas/patient_cases.xlsx")
for case_nr=1:3
    pdiag = file["Sheet$case_nr"][1:end, 1:end]
end
close(file)
# Comparison with Alisa's pdiag in hp 
maximum(abs,posterior.-pdiag)
hcat(pdiag,posterior)

