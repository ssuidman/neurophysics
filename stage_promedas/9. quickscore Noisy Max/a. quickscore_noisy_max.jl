include("../0. useful/packages.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../0. useful/useful_functions.jl")
include("../0. useful/quickscore_algorithm.jl")

function loop_term_fast(sign::Int64,pfminneg::Vector{T},previn_pfminneg::Matrix{T},one_min_previn::Matrix{T},prod_pfmin::Matrix{T}) where {T <: Real}
    A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
    B = prod_pfmin .* pfminneg'
    term = sign .* prod(A) .* [1 ; B' ./ A']; # --> prod(A + (B-A) * I0), with 'I0' the diagonal. 
    return term
end

function quickscore_noisy_max(patientprev::Matrix{Float64}, sens_normal::Matrix{Float64}, sens_medium::Array{Float64,3}, sens_extreme::Matrix{Float64}, method::String; threading::Bool=false)
# function quickscore_noisy_max(previn::Vector{Float64}, pfmin::Matrix{Float64}, pfminneg::Vector{Float64}, method::String; threading::Bool=true)
    previn = patientprev[:,1]
    pfmin = sens_extreme'
    pfmin1 = sens_medium[:,:,1]'
    pfmin2 = sens_medium[:,:,2]'
    pfminneg = sens_normal[:,1]
    
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    m1 = size(pfmin1,1) 
    println("m = $(m+m1) (postive tests)")
    
    if      occursin("Fl32",method);    previn, pfmin, pfmin1, pfmin2, pfminneg = [Float32.(x)  for x in [previn, pfmin, pfmin1, pfmin2, pfminneg]]; 
    elseif  occursin("Fl128",method);   previn, pfmin, pfmin1, pfmin2, pfminneg = [Float128.(x) for x in [previn, pfmin, pfmin1, pfmin2, pfminneg]]; 
    elseif  occursin("MF3",method);     previn, pfmin, pfmin1, pfmin2, pfminneg = [Float64x3.(x) for x in [previn, pfmin, pfmin1, pfmin2, pfminneg]]; 
    elseif  occursin("MF4",method);     previn, pfmin, pfmin1, pfmin2, pfminneg = [Float64x4.(x) for x in [previn, pfmin, pfmin1, pfmin2, pfminneg]]; 
    elseif  occursin("BF",method);      previn, pfmin, pfmin1, pfmin2, pfminneg = [BigFloat.(x) for x in [previn, pfmin, pfmin1, pfmin2, pfminneg]]; 
    end
    float_type = typeof(previn[1])

    if occursin("thread",method); threading=true; end;
    
    pfplus = zeros(float_type,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    pfplus_matrix = zeros(float_type,2^m,2^m1,n+1,1)
    previn_pfminneg = previn' .* pfminneg' 
    one_min_previn = 1 .- previn'
    myset_matrix = collect(powerset([1:m...]))
    myset_matrix1 = collect(powerset([1:m1...]))
    power_m1 = 2^m1

    println("Method: '$method'")
    dt = @elapsed begin
        if threading 
            println("Multi-threading: $(Threads.nthreads()) threads");
            pfplus_matrix = [zeros(float_type,n+1,1) for i=1:Threads.nthreads()] 
            Threads.@threads for i in ProgressBar(1:2^m) # In VS Code --> settings (left below) --> choose settings --> type in "Julia: Num Threads" and click "Edit in settings.json" --> set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
                for j in 1:2^m1
                    prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-j+1],:],dims=1).*prod(pfmin2[myset_matrix1[j],:],dims=1)
                    sign = (-1)^(length(myset_matrix[i])+length(myset_matrix1[power_m1-j+1]))
                    pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
                end
            end 
            pfplus = sum(cat(pfplus_matrix...,dims=3),dims=3)
        else
            for j in ProgressBar(1:2^m1) # iterate over 2^m possibilities 
                for i in 1:2^m
                    prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-j+1],:],dims=1).*prod(pfmin2[myset_matrix1[j],:],dims=1)
                    sign = (-1)^(length(myset_matrix[i])+length(myset_matrix1[power_m1-j+1]))
                    pfplus_matrix[i,j,:,:] = loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
                    # pfplus_matrix[i,:,:] = loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
                    println(vcat(digits(Int(sum(2 .^myset_matrix1[j]./2)),base=2,pad=4),digits(Int(sum(2 .^myset_matrix[i]./2)),base=2,pad=3)))
                    println(i,' ',j)
                    println(Float64.((prod(pfmin1[myset_matrix1[power_m1-j+1],:],dims=1).*prod(pfmin2[myset_matrix1[j],:],dims=1))[1:5]'))
                    println(Float64.((loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin).*[1;previn])[1:6]),'\n')
                end 
            end    
            pfplus = sum(pfplus_matrix,dims=(1,2))[1,1,:,:]
        end
    end
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus_matrix, pfplus, P_joint, posterior, dt 
end 

# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(13,n_disease=1000);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_noisy_max(previn, pfmin, pfminneg,"trick Fl128 thread");  
# pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread");  



