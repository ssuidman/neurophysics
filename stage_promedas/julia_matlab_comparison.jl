# THESE ARE SOME PROFILING TOOLS
# @pprof pfplus, posterior, dt = quickscore(previn,pfmin,pfminneg,true)
# @code_warntype quickscore(previn, pfmin, pfminneg, false)
# @benchmark log.(rand(101,100))
# julia --trace-compile=stderr main.jl 


module test_module
    include("packages.jl")
    include("quickscore_preparation.jl")
    include("quickscore_algorithm.jl")

    using LoopVectorization
    
    # THIS IS A COMPARISON BETWEEN MOST IMPORTANT LINES 
    function f(method)
        Random.seed!(1234) 
        x, myset, pfmin, prevminneg, prev, M, n = rand(101, 100), [1,4,5,6], rand(18,100), rand(101,100), rand(101,100), 20000, 100
        y, z, m, z_int, c = zeros(n), rand(M,n), rand(1,n), rand([-1 1],M,n), Matrix{Float64}(undef,M,1)
        zzz = zeros(101,1)

        t = @elapsed begin
            if method==1;       println("log.(x)");                             for k=1:2^20;   y = log.(x); end;
            elseif method==1.1; println("log.(x) @avx");                        for k=1:2^20;   y = @avx log.(x); end;
            elseif method==2;   println("prod(x,dims=1)");                      for k=1:2^24;   y = prod(x,dims=1); end;
            elseif method==2.1; println("prod(x,dims=1) @avx");                 for k=1:2^24;   y = @avx prod(x,dims=1); end;
            elseif method==3;   println("QS +/- exp-sum-log");                  for k=1:2^20;   y = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2)); end;
            elseif method==3.1; println("QS +/- exp-sum-log @avx log");         for k=1:2^20;   y = ((-1)^length(myset)) .* exp.(sum(@avx(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))), dims=2)); end; 
            elseif method==3.2; println("QS +/- exp-sum-log @avx log 2 steps"); for k=1:2^20;   a = 1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)); y = ((-1)^length(myset)) .* exp.(sum(@avx(log.(a)), dims=2)); end; 
            elseif method==3.3; println("QS +/- exp-sum-log @avx small prod");  for k=1:2^17;   y = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (@avx(prod(pfmin[myset, :],dims=1)) .* prevminneg .+ (1 .- prev))), dims=2)); end;
            elseif method==3.4; println("QS +/- exp-sum-log @avx exp");         for k=1:2^17;   y = ((-1)^length(myset)) .* @avx(exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))); end;
            elseif method==4;   println("QS +/- prod");                         for k=1:2^19;   y = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); end;
            elseif method==5;   println("QS log");                              for k=1:2^17;   y = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))); end;
            elseif method==6;   println("QS no prod/log");                      for k=1:2^19;   y = (1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))); end;
            elseif method==7;   println("QS prod");                             for k=1:2^19;   y = prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); end;

            elseif method==8;   println("Vector: $(typeof(z))");                for k=1:1e3;    c = z*m'; end; 
            elseif method==9;   println("Vector: $(typeof(z_int))");            for k=1:1e3;    c = z_int*m'; end; 
            elseif method==10;  println("Vector Octavian: $(typeof(z))");       for k=1:1e3;    Octavian.matmul!(c, z, m'); end; 
            elseif method==11;  println("Vector Octavian: $(typeof(z_int))");   for k=1:1e3;    Octavian.matmul!(c, z_int, m'); end; 
            elseif method==12;  println("Loop: $(typeof(z))");                  for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
            elseif method==13;  println("Loop: $(typeof(z_int))");              for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z_int[i,j]*m[j]; end; end; end; 
            elseif method==14;  println("Loop @inbounds: $(typeof(z))");        for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z[i,j]*m[j]; end; end; end; 
            elseif method==15;  println("Loop @inbounds: $(typeof(z_int))");    for k=1:1e3;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z_int[i,j]*m[j]; end; end; end; 
            elseif method==16;  println("Loop @avx: $(typeof(z))");             for k=1:1e3;    c = zeros(1,M);    @avx for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
            end
        end
        println("$t sec")
    end

    # f(1)
    # f(1.1)
    f(2)
    f(2.1)
    # f(3)
    # f(3.1)
    # f(3.2)
    
end


# 1.    For taking log of matrix(101,100)   [test 2^20  Julia.ipynb: 42.53 sec      Matlab: 24.87 sec     Julia.jl: 41.58   Julia-@avx: 24.52]
# 2.    For taking prod of matrix(101,100)  [test 2^24  Julia.ipynb: 24.23 sec      Matlab: 74.43 sec     Julia.jl: 24.71   Julia-@avx: 24.50] 
# 3.    Quickscore total log                [test 2^20  Julia.ipynb: 62.19 sec      Matlab: 39.98 sec     Julia.jl: 61.78   Julia-@avx: 31.49   Julia-@avx-2steps: 39.54]
# 4.    Quickscore total prod               [test 2^22  Julia.ipynb: 49.78 sec      Matlab: 47.48 sec     Julia.jl: 44.88]
# 5.    Quickscore including taking the log [test 2^20  Julia.ipynb: 58.45 sec      Matlab: 33.69 sec     Julia.jl: 62.51]
# 6.    Quickscore excluding taking the log [test 2^22  Julia.ipynb: 37.33 sec      Matlab: 35.13 sec     Julia.jl: 39.13]
# 7.    Quickscore taking the prod          [test 2^22  Julia.ipynb: 47.17 sec      Matlab: 47.11 sec     Julia.jl: 48.18]

#####################################################################################################################################
################################    DIFFERENT ANALYSIS: Running function in module multiple times    ################################
#####################################################################################################################################
# 8.    Matrix multiplication vector            [test 1e30       terminal:  3.63 sec,  3.54 sec   3.55 sec,  3.50 sec,  3.55 sec,  3.52 sec]
# 10.   Matrix multiplication loop              [test 1e30       terminal: 19.65 sec, 19.69 sec, 19.58 sec, 20.16 sec, 19.59 sec, 19.59 sec]
# 11.   Matrix multiplication loop @inbounds    [test 1e30       terminal: 15.48 sec, 15.54 sec, 15.69 sec, 15.74 sec, 15.71 sec, 15.50 sec]


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

# 1. No loop:   0.60
# 2. No loop:   0.57
# 3. No loop:   0.56
# 4. No loop:   0.56
# 5. No loop:   0.55
# 6. No loop:   0.56
# 7. Loop:      0.91 --> without @inbound
# 8. Loop:      0.91 --> without @inbound
# 9. Loop:      0.91 --> without @inbound
# 10.Loop:      0.92 --> without @inbound
# 11.Loop:      0.90 --> without @inbound
# 12.Loop:      0.90 --> without @inbound



