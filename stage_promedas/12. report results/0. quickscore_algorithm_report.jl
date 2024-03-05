include("../0. useful/packages.jl")

function loop_term(method::String,myset::Vector{Int64},previn::Vector{T},pfmin::Matrix{T},pfminneg::Vector{T},prevminneg::Matrix{T},prev::Matrix{T},previn_pfminneg::Matrix{T},one_min_previn::Matrix{T},pfminneg_min_previn::Matrix{T},prod_pfmin::Matrix{T}) where {T <: Real}
    if occursin("prod",method);             term = ((-1)^length(myset)) .* prod(prod_pfmin .* prevminneg .+ (1 .- prev),dims=2); 
    elseif occursin("exp-sum-log",method);  term = ((-1)^length(myset)) .* exp.(sum(log.(prod_pfmin .* prevminneg .+ (1 .- prev)), dims=2));      
    # elseif  occursin("prod",method) && occursin("entire term",method);      term = BigFloat.(((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ prev_min_1,dims=2));
    elseif occursin("trick",method);      
        A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
        B = prod_pfmin .* pfminneg'
        term = ((-1)^length(myset)) .* prod(A) .* [1 ; B' ./ A']; # --> prod(A + (B-A) * I0), with I0=[zeros(n) I(n)]' the diagonal. 
    end
    return term
end

function quickscore_report(previn::Vector{Float64}, pfmin::Matrix{Float64}, pfminneg::Vector{Float64}, method::String; threading::Bool=false)
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $m (postive tests)")
    
    if      occursin("Fl32",method);                            previn, pfmin, pfminneg = [Float32.(x)  for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("Fl128",method) && !occursin("BF",method); previn, pfmin, pfminneg = [Float128.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("MF2",method);                             previn, pfmin, pfminneg = [Float64x2.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("MF3",method);                             previn, pfmin, pfminneg = [Float64x3.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("MF4",method);                             previn, pfmin, pfminneg = [Float64x4.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("MF5",method);                             previn, pfmin, pfminneg = [Float64x5.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("D64",method);                             previn, pfmin, pfminneg = [Double64.(x) for x in [previn, pfmin, pfminneg]]; 
    elseif  occursin("BF",method);                              previn, pfmin, pfminneg = [BigFloat.(x) for x in [previn, pfmin, pfminneg]]; 
        if  occursin("Fl128",method); setprecision(BigFloat,113); elseif occursin("Fl64",method); setprecision(BigFloat,53); else; setprecision(BigFloat,300); end; 
    end
    float_type = typeof(previn[1])
    
    pfplus = zeros(float_type,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    if m<=18
        pfplus_matrix = zeros(float_type,2^m,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    else 
        println("m>18, some methods not possible")
        pfplus_matrix = zeros(m)
    end

    if occursin("thread",method); threading = true;  end; 

    prev = previn' .+ [zeros(float_type,n) I(n)]' .* (1 .- previn')
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(float_type,n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)

    previn_pfminneg = previn' .* pfminneg' 
    one_min_previn = 1 .- previn'
    pfminneg_min_previn = pfminneg' .- previn'
    
    if m>21; println("Creating powerset..."); end; 
    if occursin("linear",method); myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)]; else; myset_matrix = [powerset(1:m)...]; end;
    
    println("Method: '$method'")
    dt_array = zeros(m+1)
    dt = @elapsed begin
        if method == "MATLAB"
            println("'This is exactly the same as the code in MATLAB")
            t1 = time(); 
            for i in 0:(2^m-1) # iterate over 2^m possibilities 
                v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
                myset = findall(v.==1) # find places of the 1-elements
                if length(myset)==0 
                    # pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ prevminneg .+ (1 .- prev)), dims=2))
                    pfplus .+= ((-1)^length(myset)) .* prod(1e-50 .+ prevminneg .+ (1 .- prev), dims=2)
                elseif length(myset)==1 
                    # pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev)), dims=2))
                    pfplus .+= ((-1)^length(myset)) .* prod(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev), dims=2)
                    t2 = time()
                    dt = t2-t1 
                    dt_array[myset[1]] = dt
                    println(myset[1],"\t",v,"\t",dt)
                else 
                    # pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
                    pfplus .+= ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2)
                end
            end
        elseif threading
            println("Multi-threading: $(Threads.nthreads()) threads");
            pfplus_matrix = [zeros(float_type,n+1,1) for i=1:Threads.nthreads()] 
            Threads.@threads for i in ProgressBar(1:2^m) # In VS Code --> settings (left below) --> choose settings --> type in "Julia: Num Threads" and click "Edit in settings.json" --> set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term(method,myset_matrix[i],previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
                # pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term(method,myset,previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
                # pfplus_matrix[i,:,:] = loop_term(method,myset,previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
            end
            pfplus_matrix = permutedims(cat(pfplus_matrix...,dims=3),(3,1,2))
        elseif occursin("dict",method)
            prod_pfmin_dict = Dict{Vector{Int64},Matrix{float_type}}(myset_matrix[i] => prod(pfmin[myset_matrix[i],:],dims=1) for i in 1:m+1); 
            for i in ProgressBar(1:2^m) # iterate over 2^m possibilities 
                if occursin("dict",method) && length(myset_matrix[i])>1
                    prod_pfmin = prod_pfmin_dict[myset_matrix[i][1:end-1]] .* prod_pfmin_dict[myset_matrix[i][[end]]]  
                    push!(prod_pfmin_dict,myset_matrix[i] => prod_pfmin)
                else
                    prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                end
                pfplus_matrix[i,:,:] = loop_term(method,myset_matrix[i],previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
            end
        elseif occursin("vector",method)
            println("This method is slower than the loop version (but faster than normal) and in addition you cannot keep track of progress or use threading...")
            v_matrix = permutedims(hcat(digits.(0:2^m-1,base=2,pad=m)...))
            prod_pfmin_matrix = prod(reshape(v_matrix,size(v_matrix)...,1) .* reshape(pfmin,1,size(pfmin)...) .+ (1 .- v_matrix), dims=2)[:,1,:]
            A = prod_pfmin_matrix .* pfminneg' .* previn' .+ (1 .- previn')
            B = (prod_pfmin_matrix .* pfminneg' .- 1) .* (1 .- previn')
            pfplus_matrix = (-1).^sum(v_matrix,dims=2) .* prod(A,dims=2) .* [1 ; (A .+ B)' ./ A']'
            pfplus_matrix = reshape(pfplus_matrix,size(pfplus_matrix)...,1)
        elseif occursin("@inbounds",method)
            for i in ProgressBar(1:2^m) # iterate over 2^m possibilities 
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                @inbounds pfplus_matrix[i,:,:] = loop_term(method,myset_matrix[i],previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
            end
        elseif occursin("matrix",method)
            for i in ProgressBar(1:2^m) # iterate over 2^m possibilities 
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                pfplus_matrix[i,:,:] = loop_term(method,myset_matrix[i],previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
            end
        else
            for i in ProgressBar(1:2^m) # iterate over 2^m possibilities 
                prod_pfmin = prod(pfmin[myset_matrix[i],:],dims=1)
                pfplus = pfplus .+ loop_term(method,myset_matrix[i],previn,pfmin,pfminneg,prevminneg,prev,previn_pfminneg,one_min_previn,pfminneg_min_previn,prod_pfmin)
            end
        end
    end
    dt_array[end] = dt
    dt_array = dt_array[2:end]
    pfplus = sum(pfplus_matrix)!=0 ? sum(pfplus_matrix,dims=1)[1,:,:] : pfplus
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus_matrix, pfplus, P_joint, posterior, dt, dt_array
end 
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_report(previn, pfmin, pfminneg,"trick Fl128");  

