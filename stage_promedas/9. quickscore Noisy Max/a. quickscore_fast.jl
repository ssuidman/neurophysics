include("../0. useful/packages.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../0. useful/useful_functions.jl")
include("../0. useful/quickscore_algorithm.jl")

function loop_term_fast(myset::Vector{Int64},pfminneg::Vector{T},previn_pfminneg::Matrix{T},one_min_previn::Matrix{T},prod_pfmin::Matrix{T}) where {T <: Real}
    A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
    B = prod_pfmin .* pfminneg'
    term = ((-1)^length(myset)) .* prod(A) .* [1 ; B' ./ A']; # --> prod(A + (B-A) * I0), with 'I0' the diagonal. 
    return term
end

function quickscore_fast(previn::Vector{Float64}, pfmin::Matrix{Float64}, pfminneg::Vector{Float64}, method::String; threading::Bool=true)
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $m (postive tests)")
    
    if      occursin("Fl32",method);    previn, pfmin, pfminneg = [Float32.(x)  for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("Fl128",method);   previn, pfmin, pfminneg = [Float128.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("BF",method);      previn, pfmin, pfminneg = [BigFloat.(x) for x in [previn, pfmin, pfminneg]]; 
    end
    float_type = typeof(previn[1])
    
    pfplus = zeros(float_type,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    pfplus_matrix = zeros(float_type,2^m,n+1,1)
    previn_pfminneg = previn' .* pfminneg' 
    one_min_previn = 1 .- previn'
    myset_matrix = collect(powerset([1:m...]))

    println("Method: '$method'")
    dt = @elapsed begin
        if threading 
            println("Multi-threading: $(Threads.nthreads()) threads");
            pfplus_matrix = [zeros(float_type,n+1,1) for i=1:Threads.nthreads()] 
            Threads.@threads for i in ProgressBar(1:2^m) # In VS Code --> settings (left below) --> choose settings --> type in "Julia: Num Threads" and click "Edit in settings.json" --> set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term_fast(myset_matrix[i],pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
            end
            pfplus = sum(cat(pfplus_matrix...,dims=3),dims=3)
        else
            for i in ProgressBar(1:2^m) # iterate over 2^m possibilities 
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                pfplus_matrix[i,:,:] = loop_term_fast(myset_matrix[i],pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
            end    
            pfplus = sum(pfplus_matrix,dims=1)[1,:,:]
        end
    end
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus_matrix, pfplus, P_joint, posterior, dt 
end 

previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(13,n_disease=1000);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_fast(previn, pfmin, pfminneg,"trick Fl128 thread");  
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread");  
maximum(abs,posterior .- posterior_)



