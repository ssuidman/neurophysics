include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../4. copying MATLAB code/b. copying_functions.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")
include("../9. quickscore Noisy Max/a. quickscore_noisy_max.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"MATLAB"); 
k = 7
cases = ["2","3","6","9","42","125","128","161","404_case_1","404_case_2","404_case_3"]
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"; 
matfile = matopen(joinpath(matlab_dir,"variables/noisy_max/patient_$(cases[k]).mat"),"r"); 
    patientprev = read(matfile, "patientprev"); 
    sens_normal = read(matfile, "sens_normal"); 
    sens_medium = read(matfile, "sens_medium"); 
    sens_extreme = read(matfile, "sens_extreme"); 
    diagset = Int.(read(matfile, "diagset")); 
    pdiag_sam = read(matfile, "pdiag"); 
close(matfile)

# pfmin = vcat(prod(sens_medium,dims=3)[:,:,1]',sens_extreme')
# pfminneg = copy(sens_normal)
# previn = copy(patientprev)
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_noisy_max(patientprev,sens_normal,sens_medium,sens_extreme,"trick BF");
# (pfplus_matrix[2,15,:,1]'.*[1;previn])[1:6]'
maximum(abs,posterior.-pdiag_sam)

file = XLSX.openxlsx("/Users/sam/Large_files/Data/stage_promedas/patient_case_eabove.xlsx")
    if occursin("128",cases[k])
        pdiag = Float64.(file["Sheet4"][7:226, 5])[:,1]
        P_joint_alisa = Float64.(file["Sheet4"][7:226, 3])[:,1]
    elseif occursin("125",cases[k])
        pdiag = Float64.(file["Sheet5"][7:405, 5])[:,1]
        P_joint_alisa = Float64.(file["Sheet5"][7:405, 3])[:,1]
    end
close(file)

# # patient 404 
if occursin("case",cases[k])
    patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case "*cases[k][end]);
    pdiag = data_alisa[:,"Pdiag (high precision)"]
    P_joint_alisa = data_alisa[:,"Pjoint (high precision)"]
end

# Comparison with Alisa's pdiag in hp 
maximum(abs,P_joint.-P_joint_alisa)
maximum(abs,posterior.-pdiag)
hcat(P_joint,P_joint_alisa)
hcat(posterior,pdiag)

df = DataFrame(
    "diagset"=>diagset[:,1],
    "pdiag"=>posterior,
    "P_joint"=>P_joint,
    "P_findings"=>pfplus[1,1],
    "patientprev_1"=>patientprev[:,1],
    "sens_normal_1"=>sens_normal[:,1],
    "sens_medium_1_1"=>sens_medium[:,1,1],
    "sens_medium_1_2"=>sens_medium[:,1,2],
    "sens_medium_2_1"=>sens_medium[:,2,1],
    "sens_medium_2_2"=>sens_medium[:,2,2],
    "sens_medium_3_1"=>sens_medium[:,3,1],
    "sens_medium_3_2"=>sens_medium[:,3,2],
    "sens_medium_4_1"=>sens_medium[:,4,1],
    "sens_medium_4_2"=>sens_medium[:,4,2],
)

XLSX.writetable("/Users/sam/Downloads/patient_128_output_sam_3.xlsx",df)
CSV.write("/Users/sam/Downloads/patient_128_output_sam_3.csv",df)

# m = 6
# myset_matrix = collect(powerset([1:m...]))
# for i=1:size(myset_matrix,1)
#     println(myset_matrix[2^m-i+1],myset_matrix[i])
#     # println(sort([reverse(myset_matrix)[i];myset_matrix[i]]))
# end



