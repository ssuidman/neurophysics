include("../0. useful/packages.jl")
include("b. time_comparison_function.jl")

println("Total running time will take ~3.5 minutes for all methods + ~10.5 minutes for 'log.(x_matrix) BF Fl128'") 
# Change 'method_names' and 'method_numbers' in Julia and 'methods' in matlab when adding a new line in 'b. time_comparison_function.jl'
method_names = ["log.(x_matrix)","@avx(log.(x_matrix))","log.(x_matrix) BF Fl128 (n=10^11)","prod(x_matrix,dims=1)","@avx(prod(x_matrix,dims=1))","QS exp-sum-log","QS exp-sum-log @avx(log)","QS exp-sum-log @avx(log) 2 steps","QS exp-sum-log @avx(prod(pfmin))","QS exp-sum-log @avx(exp)","QS exp-sum-log broadcast","QS exp-sum-log declare-type","QS prod","QS prod @avx(prod(term))","QS log","QS log @avx","QS no prod/log","QS no prod/log @avx","Vector: Matrix{Float64}","Vector: Matrix{Int64}","Vector Octavian: Matrix{Float64}","Vector Octavian: Matrix{Int64}","Loop: Matrix{Float64}","Loop: Matrix{Int64}","Loop @inbounds: Matrix{Float64}","Loop @inbounds: Matrix{Int64}","Loop @avx: Matrix{Float64}","Loop @avx: Matrix{Int64}"];
method_numbers = [1,1.1,1.2,2,2.1,3,3.1,3.2,3.3,3.4,3.5,3.6,4,4.1,5,5.1,6,6.1,7,7.1,8,8.1,9,9.1,10,10.1,11,11.1];

dt_julia = zeros(size(method_numbers))
for (i,method) in enumerate(method_numbers)
    dt_julia[i] = f(method) 
end 

# run(`/Applications/MATLAB_R2022b.app/bin/matlab -nodesktop -nosplash -r "run('~/Documents/MATLAB/Wim/2. running time comparison/c_comparison_run.m');exit;"`)
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables/running_time_comparison/matlab_times.mat"),"r")
    dt_matlab = read(matfile, "dt_results")[1,:];
close(matfile)

df = DataFrame([method_names,dt_julia,dt_matlab],["Methods","Julia","MATLAB"])
CSV.write("variables/time_comparison/time_method_comparison.csv",df)


# 1.    log of matrix(101,100)              [test 2^20  Julia.ipynb: 42.53 sec      Matlab: 24.87 sec     Julia.jl: 41.58   Julia-@avx: 24.52]
# 2.    prod of matrix(101,100)             [test 2^24  Julia.ipynb: 24.23 sec      Matlab: 74.43 sec     Julia.jl: 24.71   Julia-@avx: 24.50] 
# 3.    Quickscore +/- exp-sum-log.()       [test 2^20  Julia.ipynb: 62.19 sec      Matlab: 39.98 sec     Julia.jl: 61.78   Julia-@avx: 31.49   Julia-@avx-2steps: 39.54]
# 4.    Quickscore +/- prod()               [test 2^22  Julia.ipynb: 49.78 sec      Matlab: 47.48 sec     Julia.jl: 44.88   Julia-@avx: 43.76]
# 5.    Quickscore log.()                   [test 2^20  Julia.ipynb: 58.45 sec      Matlab: 33.69 sec     Julia.jl: 62.51   Julia-@avx: 27.74]
# 6.    Quickscore ()                       [test 2^22  Julia.ipynb: 37.33 sec      Matlab: 35.13 sec     Julia.jl: 39.13   Julia-@avx: 34.80]
# 7.    Quickscore prod()                   [test 2^22  Julia.ipynb: 47.17 sec      Matlab: 47.11 sec     Julia.jl: 48.18   Julia-@avx: 42.82]

#####################################################################################################################################
################################    DIFFERENT ANALYSIS: Running function in module multiple times    ################################
#####################################################################################################################################
# 8.    Matrix multiplication vector            [test 1e30       terminal:  3.63 sec,  3.54 sec   3.55 sec,  3.50 sec,  3.55 sec,  3.52 sec]
# 10.   Matrix multiplication loop              [test 1e30       terminal: 19.65 sec, 19.69 sec, 19.58 sec, 20.16 sec, 19.59 sec, 19.59 sec]
# 11.   Matrix multiplication loop @inbounds    [test 1e30       terminal: 15.48 sec, 15.54 sec, 15.69 sec, 15.74 sec, 15.71 sec, 15.50 sec]

######## Even een bevinding 
# Nested for loops zijn heel erg langzaam in VS Code als het niet een functie zit. Anders is het ongeveer even snel volgens mij. 
# Het maakt hierbij niet uit of het met @inbounds is of niet. 
# Dit ook kijken in Matlab 
# Ook kijken naar andere dingen die ik tot nu toe heb gedaan in Matlab 
# Beginnen om ergens een overzicht te maken 


###############################################################################################################################
###################################################   Results test_giel.jl  ###################################################
###############################################################################################################################
###  Remarks: There is no precompilation during the minute of running time via "julia --trace-compile=stderr test_giel.jl"  ###
###############################################################################################################################
# 1. Loop:      61.49     5.74    60.48
# 2. No loop:   68.86     6.50    67.78
# 3. Loop:      59.62     5.76    58.29
# 4. Loop:      47.95     5.80    47.55
# 5. Loop:      47.83     5.80    48.00
# 6. No loop:   55.11     6.56    55.45
# 7. No loop:   54.12     6.46    54.85
# 8. No loop:   55.29     6.57    55.52

# 1. Loop:      61.95   5.69  60.90
# 2. No loop:   69.11   6.32  69.04
# 3. Loop:      59.88   5.58  58.99
# 4. No loop:   55.16   6.33  54.43
# 5. Loop:      47.82   5.60  47.71
# 6. No loop:   55.55   6.40  55.27
# 7. Loop:      47.93   5.74  47.75
# 8. No loop:   55.79   6.44  55.66

# 1. Loop:   61.05   61.70
# 2. Loop:   61.50   62.14
# 3. Loop:   58.58   59.71
# 4. Loop:   47.63   47.68
# 5. Loop:           47.78
# 6. Loop:           48.12

# 1. No loop:    66.76  68.78
# 2. No loop:    67.97  68.84
# 3. No loop:    65.38  65.80
# 4. No loop:    54.37  54.48
# 5. No loop:    55.36  55.06
# 6. No loop:    55.76  55.46

# 1. Loop:      61.08
# 2. Loop:      61.87
# 3. Loop:      58.96
# 4. No loop:   54.48
# 5. No loop:   54.18
# 6. No loop:   55.18

# 1. No loop:   67.40
# 2. No loop:   68.89
# 3. No loop:   66.02
# 4. Loop:      47.69
# 5. Loop:      47.95
# 6. Loop:      48.10

######## with Module #######
# 1. Loop @inbounds     0.50    5.40
# 2. Loop @inbounds     0.49    5.43
# 3. Loop @inbounds     0.50    5.47
# 4. Loop @inbounds     0.47    5.59
# 5. Loop @inbounds     0.47    5.60
# 6. Loop @inbounds     0.47    5.59
# 7. Vector built-in    0.54    6.35
# 8. Vector built-in    0.53    6.28
# 9. Vector built-in    0.55    6.24
# 10. Vector built-in   0.55    6.35
# 11. Vector built-in   0.55    6.44
# 12. Vector built-in   0.54    6.45

######## with Module #######
# 1. Loop no @inbounds  0.99    10.00
# 2. Loop no @inbounds  0.95    9.87
# 3. Loop no @inbounds  0.95    10.07
# 4. Loop no @inbounds  0.93    10.09
# 5. Loop no @inbounds  0.93    10.19
# 6. Loop no @inbounds  0.95    10.14
# 7. Vector Octavian    0.41    5.05
# 8. Vector Octavian    0.42    4.87
# 9. Vector Octavian    0.42    5.02
# 10. Vector Octavian   0.42    4.96
# 11. Vector Octavian   0.42    4.93
# 12. Vector Octavian   0.42    4.92 

###### without Module #####
# 1. Loop no @inbounds  0.51    5.50
# 2. Loop no @inbounds  0.47    5.61
# 3. Loop no @inbounds  0.49    5.63
# 4. Loop no @inbounds  0.47    5.57
# 5. Loop no @inbounds  0.47    5.71
# 6. Loop no @inbounds  0.48    5.67
# 7. Vector built-in    0.56    6.30
# 8. Vector built-int   0.55    6.37
# 9. Vector built-int   0.55    6.32
# 10. Vector built-int  0.56    6.34
# 11. Vector built-int  0.55    6.42
# 12. Vector built-int  0.55    6.38

###### without Module ##### 
# 1. Loop no @inbounds 1.00     10.06
# 2. Loop no @inbounds 0.96     10.13
# 3. Loop no @inbounds 0.96     10.07
# 4. Loop no @inbounds 0.94     10.11
# 5. Loop no @inbounds 0.94     10.07
# 6. Loop no @inbounds 0.97     10.07
# 7. Vector Octavian 0.43       4.98
# 8. Vector Octavian 0.44       5.01
# 9. Vector Octavian 0.43       4.94
# 10. Vector Octavian 0.44      5.00
# 11. Vector Octavian 0.43      5.04
# 12. Vector Octavian 0.43      4.96




