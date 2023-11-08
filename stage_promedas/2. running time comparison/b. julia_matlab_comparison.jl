# THESE ARE SOME PROFILING TOOLS
# @pprof pfplus, posterior, dt = quickscore(previn,pfmin,pfminneg,true)
# @code_warntype quickscore(previn, pfmin, pfminneg, false)
# @benchmark log.(rand(101,100))
# julia --trace-compile=stderr main.jl 

include("../packages.jl") 
include("../quickscore_algorithm.jl") 

# THIS IS A COMPARISON BETWEEN MOST IMPORTANT LINES 
function f(method) 
    """
    This function aims to compare different one-line-methods in order to compare different runnings times after running the 
    test_giel file. This could be beneficial for the actual quickscore algorithm in the end. Some conclusions are 
    1)  'prod' is much faster than 'exp-sum-log'
    2)  @avx can speed up things, but it can give another precision (therefore @avx(log) 2 steps is needed), so you need to be extremely
        carefull when using this. Other methods that can speed up things are @inbounds (used to remove checks on loop bounds). 
    3) Declaring types in the function (3.6) is not fast if you do it in line. 
    4) Declaring variables can help in the precompiling, but is not important for the loop part at all. 
    """
    Random.seed!(1234) 
    x_matrix, myset, pfmin, prevminneg, prev, M, n = rand(101, 100), [1,4,5,6], rand(18,100), rand(101,100), rand(101,100), 200000, 100 
    y, z, m, z_int, c = zeros(n), rand(M,n), rand(1,n), rand([-1 1],M,n), Matrix{Float64}(undef,M,1) 

    t = @elapsed begin 
        if method==1;       println("log.(x_matrix)");                      for k=1:2^18;   y = log.(x_matrix); end; 
        elseif method==1.1; println("@avx(log.(x_matrix))");                for k=1:2^18;   y = @avx(log.(x_matrix)); end; 
        elseif method==2;   println("prod(x_matrix,dims=1)");               for k=1:2^18;   y = prod(x_matrix,dims=1); end; 
        elseif method==2.1; println("@avx(prod(x_matrix,dims=1))");         for k=1:2^18;   y = @avx(prod(x_matrix,dims=1)); end; 
        elseif method==3;   println("QS exp-sum-log");                      for k=1:2^18;   y = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2)); end;
        elseif method==3.1; println("QS exp-sum-log @avx(log)");            for k=1:2^18;   y = ((-1)^length(myset)) .* exp.(sum(@avx(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))), dims=2)); end; 
        elseif method==3.2; println("QS exp-sum-log @avx(log) 2 steps");    for k=1:2^18;   a = 1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)); y = ((-1)^length(myset)) .* exp.(sum(@avx(log.(a)), dims=2)); end; 
        elseif method==3.3; println("QS exp-sum-log @avx(prod(pfmin))");    for k=1:2^18;   y = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (@avx(prod(pfmin[myset, :],dims=1)) .* prevminneg .+ (1 .- prev))), dims=2)); end;
        elseif method==3.4; println("QS exp-sum-log @avx(exp)");            for k=1:2^18;   y = ((-1)^length(myset)) .* @avx(exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))); end;
        elseif method==3.5; println("QS exp-sum-log broadcast");            for k=1:2^18;   y = ((-1)^length(myset)) .* exp.(sum(broadcast(log, 1e-50 .+ (prod(pfmin[myset, :], dims=1) .* prevminneg .+ (1 .- prev))), dims=2)); end; 
        elseif method==3.6; println("QS exp-sum-log declare-type");         for k=1:2^18;   y = ((-1)^length(myset::Vector{Int64})) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset::Vector{Int64}, :],dims=1) .* prevminneg::Matrix{Float64} .+ (1 .- prev))), dims=2)); end; 
        elseif method==4;   println("QS prod");                             for k=1:2^18;   y = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); end;
        elseif method==4.1; println("QS prod @avx(prod(term))");            for k=1:2^18;   y = ((-1)^length(myset)) .* @avx(prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)); end;
        elseif method==5;   println("QS log");                              for k=1:2^18;   y = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))); end;
        elseif method==5.1; println("QS log @avx");                         for k=1:2^18;   y = @avx(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))); end;
        elseif method==6;   println("QS no prod/log");                      for k=1:2^18;   y = (1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))); end;
        elseif method==6.1; println("QS no prod/log @avx");                 for k=1:2^18;   y = @avx(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))); end;
        elseif method==7;   println("Vector: $(typeof(z))");                for k=1:1e3;    c = z*m'; end; 
        elseif method==7.1; println("Vector: $(typeof(z_int))");            for k=1:1e3;    c = z_int*m'; end; 
        elseif method==8;   println("Vector Octavian: $(typeof(z))");       for k=1:1e3;    Octavian.matmul!(c, z, m'); end; 
        elseif method==8.1; println("Vector Octavian: $(typeof(z_int))");   for k=1:1e3;    Octavian.matmul!(c, z_int, m'); end; 
        elseif method==9;   println("Loop: $(typeof(z))");                  for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==9.1; println("Loop: $(typeof(z_int))");              for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z_int[i,j]*m[j]; end; end; end; 
        elseif method==10;  println("Loop @inbounds: $(typeof(z))");        for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==10.1;println("Loop @inbounds: $(typeof(z_int))");    for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z_int[i,j]*m[j]; end; end; end; 
        elseif method==11;  println("Loop @avx: $(typeof(z))");             for k=1:1e3;    c = zeros(1,M);    @avx for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==11.1;println("Loop @avx: $(typeof(z_int))");         for k=1:1e3;    c = zeros(1,M);    @avx for i=1:M; for j=1:n; c[i]+=z_int[i,j]*m[j]; end; end; end; 
        end 
    end 
    println("\t$t sec") 
end 

println("Total running time will take ~3.5 minutes") 
methods = [1,1.1,2,2.1,3,3.1,3.2,3.3,3.4,3.5,3.6,4,4.1,5,5.1,6,6.1,7,7.1,8,8.1,9,9.1,10,10.1,11,11.1]; 
for method in methods 
    f(method) 
end 



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




