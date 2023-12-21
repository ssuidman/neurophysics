include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../4. copying MATLAB code/b. copying_functions.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")
include("../9. quickscore Noisy Max/a. quickscore_noisy_max.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"MATLAB"); 
k = 8 # 128(8), 323(11)
cases = ["2","3","6","9","42","117","125","128","161","282","323","404_case_1","404_case_2","404_case_3"]; println("patient ",cases[k]);
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"; 
matfile = matopen(joinpath(matlab_dir,"variables/noisy_max/patient_$(cases[k]).mat"),"r"); 
    patientprev = read(matfile, "patientprev"); 
    sens_normal = read(matfile, "sens_normal"); 
    sens_medium = read(matfile, "sens_medium"); 
    sens_extreme = read(matfile, "sens_extreme"); 
    diagset = Int.(read(matfile, "diagset")); 
    pdiag_sam = read(matfile, "pdiag"); 
close(matfile)

file = XLSX.openxlsx("/Users/sam/Large_files/Data/stage_promedas/patient_case_eabove.xlsx")
    if occursin("128",cases[k])
        pdiag = Float64.(file["Sheet4"][5:224, 5])[:,1]
        P_joint_alisa = Float64.(file["Sheet4"][5:224, 3])[:,1]
        pdiag_higher_precision = Float64.(file["Sheet4"][5:224, 6])[:,1]
    elseif occursin("125",cases[k])
        pdiag = Float64.(file["Sheet5"][7:405, 5])[:,1]
        P_joint_alisa = Float64.(file["Sheet5"][7:405, 3])[:,1]
    elseif occursin("161",cases[k])
        pdiag = Float64.(file["Sheet6"][8:783, 4])[:,1]
        P_joint_alisa = Float64.(file["Sheet6"][8:783, 3])[:,1]
    elseif occursin("282",cases[k])
        pdiag = Float64.(file["Sheet7"][8:499, 4])[:,1]
        P_joint_alisa = Float64.(file["Sheet7"][8:499, 3])[:,1]
    elseif occursin("323",cases[k])
        pdiag = Float64.(file["Sheet8"][5:612, 4])[:,1]
        P_joint_alisa = Float64.(file["Sheet8"][5:612, 3])[:,1]
    elseif occursin("117",cases[k])
        pdiag = Float64.(file["Sheet9"][5:538, 4])[:,1]
        P_joint_alisa = Float64.(file["Sheet9"][5:538, 3])[:,1]
    end
close(file)

# # patient 404 
if occursin("case",cases[k])
    patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case "*cases[k][end]);
    pdiag = data_alisa[:,"Pdiag (high precision)"]
    P_joint_alisa = data_alisa[:,"Pjoint (high precision)"]
end

pfmin = vcat(prod(sens_medium,dims=3)[:,:,1]',sens_extreme')
pfminneg = copy(sens_normal)
previn = copy(patientprev)
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_noisy_max(patientprev,sens_normal,sens_medium,sens_extreme,"BF");

# Comparison with Alisa's pdiag in hp 
maximum(abs,P_joint.- P_joint_alisa) 
maximum(abs,posterior.- pdiag_higher_precision) 
hcat(P_joint,P_joint_alisa) 
hcat(posterior,pdiag) 
# maximum(abs,pdiag_sam.-pdiag) 

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

# XLSX.writetable("/Users/sam/Downloads/patient_128_output_sam_3.xlsx",df)
# CSV.write("/Users/sam/Downloads/patient_128_output_sam_3.csv",df)

# m = 6
# myset_matrix = collect(powerset([1:m...]))
# for i=1:size(myset_matrix,1)
#     println(myset_matrix[2^m-i+1],myset_matrix[i])
#     # println(sort([reverse(myset_matrix)[i];myset_matrix[i]]))
# end

previn, pfmin, pfmin1, pfmin2, pfminneg = patientprev[:,1], sens_extreme', sens_medium[:,:,1]', sens_medium[:,:,2]', sens_normal[:,1]; 
(m,n),m1 = (!isempty(pfmin) ? size(pfmin) : (0, length(previn))), size(pfmin1,1) 

previn_pfminneg, one_min_previn, myset_matrix, myset_matrix1, power_m1 = previn' .* pfminneg', 1 .- previn', collect(powerset([1:m...])), collect(powerset([1:m1...])), 2^m1;
i,j = 2^m, 2^m1
prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-j+1],:],dims=1).*prod(pfmin2[myset_matrix1[j],:],dims=1)
sign = (-1)^(length(myset_matrix[i])+length(myset_matrix1[power_m1-j+1]))

A = prod_pfmin .* previn_pfminneg .+ one_min_previn
B = prod_pfmin .* pfminneg'
term = sign .* prod(A) .* [1 ; B' ./ A']; # --> prod(A + (B-A) * I0), with 'I0' the diagonal. 

pfplus_matrix[i,j,:,:] = loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
# pfplus_matrix[i,:,:] = loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
println(vcat(digits(Int(sum(2 .^myset_matrix1[j]./2)),base=2,pad=4),digits(Int(sum(2 .^myset_matrix[i]./2)),base=2,pad=3)))
# println(i,' ',j)
println(Float64.((prod(pfmin1[myset_matrix1[power_m1-j+1],:],dims=1).*prod(pfmin2[myset_matrix1[j],:],dims=1))[1:5]'))
println(Float64.((loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin).*[1;previn])[1:6]),'\n')

findall(patientprev[:,2] .- one_min_previn'.!=0)
# for patient 323
idx = vcat([[i[1]-1,i[1],i[1]+1] for i in findall(patientprev[:,2] .- one_min_previn'.!=0)]...)
df = DataFrame(
    "diagn" => idx,
    "patientprev[:,1]" => patientprev[idx,1], 
    "previn" => previn[idx],
    "patientprev[:,2]" => patientprev[idx,2], 
    "one_min_previn" => one_min_previn[idx]
)





