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
    # Change the 'method_numbers' and 'method_names' in Julia and 'methods' in MATLAB by adding the new numbers you add
    # Then run MATLAB and thereafter Julia '2. running_time_comparison' again 
    t = @elapsed begin 
        if method==1;       println("log.(x_matrix)");                      for k=1:2^18;   y = log.(x_matrix); end; 
        elseif method==1.1; println("@avx(log.(x_matrix))");                for k=1:2^18;   y = @avx(log.(x_matrix)); end; 
        elseif method==1.2; println("log.(x_matrix) BF Fl128");             for k in ProgressBar(1:2^11);   x_matrix = BigFloat.(x_matrix,precision=113); y = (log.(x_matrix)); end; 
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
        elseif method==7;   println("Vector: $(typeof(z))");                for k=1:1e2;    c = z*m'; end; 
        elseif method==7.1; println("Vector: $(typeof(z_int))");            for k=1:1e2;    c = z_int*m'; end; 
        elseif method==8;   println("Vector Octavian: $(typeof(z))");       for k=1:1e2;    Octavian.matmul!(c, z, m'); end; 
        elseif method==8.1; println("Vector Octavian: $(typeof(z_int))");   for k=1:1e2;    Octavian.matmul!(c, z_int, m'); end; 
        elseif method==9;   println("Loop: $(typeof(z))");                  for k=1:1e2;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==9.1; println("Loop: $(typeof(z_int))");              for k=1:1e2;    c = zeros(1,M);    for i=1:M; for j=1:n; c[i]+=z_int[i,j]*m[j]; end; end; end; 
        elseif method==10;  println("Loop @inbounds: $(typeof(z))");        for k=1:1e2;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==10.1;println("Loop @inbounds: $(typeof(z_int))");    for k=1:1e2;    c = zeros(1,M);    for i=1:M; for j=1:n; @inbounds c[i]+=z_int[i,j]*m[j]; end; end; end; 
        elseif method==11;  println("Loop @avx: $(typeof(z))");             for k=1:1e2;    c = zeros(1,M);    @avx for i=1:M; for j=1:n; c[i]+=z[i,j]*m[j]; end; end; end; 
        elseif method==11.1;println("Loop @avx: $(typeof(z_int))");         for k=1:1e2;    c = zeros(1,M);    @avx for i=1:M; for j=1:n; c[i]+=z_int[i,j]*m[j]; end; end; end; 
        end 
    end 
    println("\t$t sec") 
    return t
end 

