# This file consists of dictionaries that contain the 10 patient cases, such that you can call using (e.g.) previn_dict["1"] for patient 1 (10 patients in total)
using MAT
using DelimitedFiles
using Quadmath
using MultiFloats
using DataFrames
using Combinatorics
using ProgressBars
using LinearAlgebra
using CSV

# This is the folder to where this file is
folder_path = "/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/11. transfer_to_alisa/"
# This is the input data from MATLAB where I created previn(prevalences), pfmin(sens_extreme for binary), pfminneg(sens_normal for binary) 
previn_dict, pfmin_dict, pfminneg_dict = Dict{String,Vector}(), Dict{String,Matrix}(), Dict{String,Vector}()
# A patient file was created in MATLAB with the following specifications 
for i=1:10
    matfile = matopen(folder_path*"extrapolation_random_patients/patient_$i.mat","r"); 
        previn_dict["$i"] = read(matfile, "previn")[1,:]; 
        pfmin_dict["$i"] = read(matfile, "pfmin"); (m,n)=size(pfmin_dict["$i"]);
        pfminneg_dict["$i"] = read(matfile, "pfminneg")[1,:]; 
        m_min = 2
        println("patient $i:\t n = $n \t m+ = $m \t m- = $m_min")
    close(matfile)
end

# BigFloat is called 'extreme', because I set the precision higher than the normal settings. pfplus BF is the probability of findings in BigFloat precision (so not the relative error)
method_names = ["pfplus BF",["trick $(type)thread" for type in ["BF extreme ","MF6 ","MF5 ","MF4 ","MF3 ","MF2 ",""]]...]
# These are the relative errors that build up to 20 positive findings for each of the precision types 
error_df_dict = Dict("$k"=>DataFrame() for k=1:10)
# These are the running times that build up to 20 positive findings for each of the precision types 
dt_df_dict = Dict("$k"=>DataFrame() for k=1:10)

for k=1:10
    dt_df = CSV.read(folder_path*"extrapolation_random_patients/patient_$(k)_dt.csv", DataFrame)    
    dt_df_dict["$k"] = dt_df[:,method_names[2:end]]
    error_df = CSV.read(folder_path*"extrapolation_random_patients/patient_$(k)_error.csv", DataFrame)    
    error_df_dict["$k"] = error_df[:,method_names]
end

# These DataFrames look for each of the 10 patient cases when precision is still just sufficient enough (relative error=δ<0.1).
# All information is based on the precious datasets dt_df/error_df and the linear extrapolation is between 8-20 positive findings. 
# The DataFrames consist of:
# method:   corresponding precision type 
# m_max:    corresponding amount of iterations --> if m_max>20 then it is based on extrapolated data 
# eps:      corresponding machine epsilon for the precision type
# pfplus:   value of p(F+,F-) --> if m_max>20 then it is based on extrapolated data 
# dt:       running time --> if m_max>20 then it is based on extrapolated data 
# error:    relative error (should be just smaller than 0.1)
# error+1:  relative error for the next iteration (should be just bigger than 0.1)
extrapolation_df_dict = Dict("$k"=>DataFrame() for k=1:10)
for k=1:10
    extrapolation_df_dict["$k"] = CSV.read(folder_path*"extrapolation_random_patients/patient_$(k)_extrapolation.csv", DataFrame)    
end
# For patient 1 the data looks like this:
dt_df_dict["1"]
error_df_dict["1"]

extrapolation_df_dict["1"]


