include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")

# function to prepare data for patient 404 case with m=9
function prepare_patient_data(case)
    patient_cases_raw = Dict{String,Matrix{Any}}()
    case_nr = parse(Int,case[end])
    if case_nr in [1,2,3]
        file = XLSX.openxlsx("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_cases.xlsx")
            for case_nr=1:3
                patient_cases_raw["case $case_nr"] = file["Sheet$case_nr"][1:end, 1:end]
            end
        close(file)
        data_alisa = DataFrame(Float64.(patient_cases_raw[case][8:end,:]),patient_cases_raw[case][6,:])
        idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
        all_sens = Matrix(data_alisa[:,idx_sens])
        pos_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'true'"
        neg_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'false'"
        sens = all_sens[:,pos_idx]
        sensneg = all_sens[:,neg_idx]; n_min = size(sensneg,2); 
        previn = data_alisa[:,"prevalence"]
        pfmin = 1 .- sens'; (m,n) = size(pfmin); 
        pfminneg = prod(1 .- sensneg', dims=1)[1,:]
        m_min = 6; 
    elseif case_nr==4
        # A patient file was created in MATLAB with the following specifications 
        # patient(999).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978];
        # patient(999).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true"];
        # patient(999).age = 40; 
        # patient(999).gender = 'F';
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_999.mat","r"); 
            previn = read(matfile, "patientprev")[:,1]; 
            pfmin = permutedims(read(matfile, "sens_extreme")); (m,n)=size(pfmin); 
            pfminneg = read(matfile, "sens_normal")[:,1]; 
            m_min = 2
        close(matfile)
    elseif case_nr==5
        # A patient file was created in MATLAB with the following specifications 
        # patient(997).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978,2687,2689,2690,2691];
        # patient(997).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true","true","true","true","true"];
        # patient(997).age = 40; 
        # patient(997).gender = 'F';
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_997.mat","r"); 
            previn = read(matfile, "previn")[1,:]; 
            pfmin = read(matfile, "pfmin"); (m,n)=size(pfmin);
            pfminneg = read(matfile, "pfminneg")[1,:]; 
            m_min = 2
        close(matfile)
    elseif case_nr==6
        # A patient file was created in MATLAB with the following specifications 
        # patient(998).test = [55,437,987,1019,1122,1420,1731,2867,3118,3235,3284,3426,3568,3812,3978,2687,2689,2690,2691,2692,2699,2702,2704,2706];
        # patient(998).tvalue = ["true","true","true","true","true","true","true","true","true","true","true","true","false","false","true","true","true","true","true","true","true","true","true","true"];
        # patient(998).age = 40; 
        # patient(998).gender = 'F';
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_998.mat","r"); 
            previn = read(matfile, "previn")[1,:]; 
            pfmin = read(matfile, "pfmin"); (m,n)=size(pfmin);
            pfminneg = read(matfile, "pfminneg")[1,:]; 
            m_min = 2
        close(matfile)
    elseif case_nr==7
        println("Noisy MAX case with m_tot=18")
        # A patient file was created in MATLAB with the following specifications 
        # patient(320) --> noisy max, but this case gives NaN because source data is not right 
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_320.mat","r"); 
            patientprev = read(matfile, "patientprev")
            sens_normal = read(matfile, "sens_normal"); 
            sens_medium = read(matfile, "sens_medium")
            sens_extreme = read(matfile, "sens_extreme")
            n, m, m_min = 1097, 18, 9
        close(matfile)
    elseif case_nr==8
        println("Noisy MAX case with m_tot=18")
        # A patient file was created in MATLAB with the following specifications 
        # patient(320) --> noisy max
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_105.mat","r"); 
            patientprev = read(matfile, "patientprev")
            sens_normal = read(matfile, "sens_normal"); 
            sens_medium = read(matfile, "sens_medium")
            sens_extreme = read(matfile, "sens_extreme")
            n, m, m_min = 695, 18, 25
        close(matfile)
    elseif case_nr==9
        println("Whole set of random patients")
        previn_dict, pfmin_dict, pfminneg_dict = Dict{String,Vector}(), Dict{String,Matrix}(), Dict{String,Vector}()
        # A patient file was created in MATLAB with the following specifications 
        for i=1:10
            matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables//extrapolation_random_patients/patient_$i.mat","r"); 
                previn_dict["$i"] = read(matfile, "previn")[1,:]; 
                pfmin_dict["$i"] = read(matfile, "pfmin"); (m,n)=size(pfmin_dict["$i"]);
                pfminneg_dict["$i"] = read(matfile, "pfminneg")[1,:]; 
                m_min = 2
                println("n = $n \t m+ = $m \t m- = $m_min")
            close(matfile)
        end
    end

    # m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    # prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    # for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    # prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    
    if case_nr in 1:6
        println("n = $n \t m+ = $m \t m- = $m_min")
        return previn, pfmin, pfminneg
    elseif case_nr in 7:8
        println("n = $n \t m+ = $m \t m- = $m_min")
        return patientprev,sens_normal,sens_medium,sens_extreme
    elseif case_nr in 9
        return previn_dict,pfmin_dict,pfminneg_dict
    end
end 
# previn, pfmin, pfminneg = prepare_patient_data("case 3"); (m,n)=size(pfmin); 
# patientprev,sens_normal,sens_medium,sens_extreme = prepare_patient_data("case 7"); 

