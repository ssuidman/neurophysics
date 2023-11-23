include("packages.jl")
include("useful_functions.jl")


function loop_term(myset::Vector{Int64},pfmin::Matrix{T},prevminneg::Matrix{T},prev::Matrix{T},method::String) where {T <: Real}
    if      method in ["exp-sum-log","exp-sum-log Fl32","exp-sum-log QM Fl128"];    term = ((-1)^length(myset)) .* exp.(sum(log.(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2));
        # elseif  method in ["exp-sum-log Fl32"];                     term = ((-1)^length(myset)) .* exp.(sum(log.(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2));
        # elseif  method in ["exp-sum-log QM Fl128"];                 term = ((-1)^length(myset)) .* exp.(sum(log.(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2));
    elseif  method in ["exp-sum-log prod1"];                    term = ((-1)^length(myset)) .* exp.(sum(log.(prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2));
    elseif  method in ["exp-sum-log BF(entire term)"];          term = BigFloat.(((-1)^length(myset)) .* exp.(sum(log.(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2)));
        # elseif  method in ["prod","prod Fl32","prod QM Fl128","prod QM Fl128 thread","prod thread"];  term = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
    elseif  occursin("prod",method);    term = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
    elseif  occursin("exp-sum-log",method);    term = ((-1)^length(myset)) .* exp.(sum(log.(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2));
    # elseif  method in ["prod Fl32"];                            term = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
    # elseif  method in ["prod QM Fl128","prod QM Fl128 thread"]; term = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
    elseif  method in ["prod prod1"];                           term = ((-1)^length(myset)) .* prod(prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2);  
    elseif  method in ["prod BF(entire term)"];                 term = BigFloat.(((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2));
    # elseif  method in ["prod BF","prod BF Fl64","prod BF Fl128","prod BF thread"];   term = ((-1)^length(myset)) .* prod(prod(BigFloat.(pfmin)[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
    # elseif  method in ["exp-sum-log BF","exp-sum-log BF thread"];               term = ((-1)^length(myset)) .* exp.(sum(log.(prod(BigFloat.(pfmin)[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2)); 
    end
    return term
end

function quickscore(previn::Vector{Float64}, pfmin::Matrix{Float64}, pfminneg::Vector{Float64}, method::String; threading::Bool=false)
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $m (postive tests)")
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    # pfplus = zeros(Float64,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    pfplus = zeros(Float64,2^m,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    
    if      occursin("Fl32",method);    previn, pfmin, prevminneg, prev, pfplus = [Float32.(x) for x in [previn, pfmin, prevminneg, prev, pfplus]]; 
    elseif  occursin("QM Fl128",method);previn, pfmin, prevminneg, prev, pfplus = [Float128.(x) for x in [previn, pfmin, prevminneg, prev, pfplus]]; 
    elseif  occursin("Fl64",method);    setprecision(BigFloat,53); 
    elseif  occursin("BF",method)       setprecision(BigFloat,256); previn, pfmin, prevminneg, prev, pfplus = [BigFloat.(x) for x in [previn, pfmin, prevminneg, prev, pfplus]]; 
    elseif  occursin("BF Fl128",method);setprecision(BigFloat,113); previn, pfmin, prevminneg, prev, pfplus = [BigFloat.(x) for x in [previn, pfmin, prevminneg, prev, pfplus]]; 
    end

    if occursin("thread",method); threading = true; end; 

    if method!="MATLAB"
        println("Method: '$method'")
        dt = @elapsed begin
            if threading
                Threads.@threads for i in ProgressBar(0:(2^m-1)) # In VS Code --> settings (left below) --> choose settings --> type in "Julia: Num Threads" and click "Edit in settings.json" --> set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
                    v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
                    myset = findall(v.==1) # find places of the 1-elements
                    pfplus[i+1,:,:] = loop_term(myset,pfmin,prevminneg,prev,method)
                    # pfplus[i+1,:,:] = rand(size(pfplus,2),size(pfplus,3))

                end
            else
                for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
                    v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
                    myset = findall(v.==1) # find places of the 1-elements
                    # pfplus = pfplus .+ loop_term(myset,pfmin,prevminneg,prev,method)
                    pfplus[i+1,:,:] = loop_term(myset,pfmin,prevminneg,prev,method)
                end
            end
        end
    else
        println("Method: 'MATLAB exp-sum-log'\t--> exactly same as Matlab code")
        dt = @elapsed begin
            for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
                v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
                myset = findall(v.==1) # find places of the 1-elements
                if length(myset)==0 
                    pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ prevminneg .+ (1 .- prev)), dims=2))
                elseif length(myset)==1 
                    pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev)), dims=2))
                else 
                    pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
                end
            end
        end
    end
    pfplus = sum(pfplus,dims=1)[1,:,:]
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus, P_joint, posterior, dt 
end 
# pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"prod");  

# previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");
# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(7,n_disease=773)


# final results
# pfplus is p(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1),
# pfpluspd is numerator in heckerman (12)  p(F+,F-),P(F+,F-|d_1=1)p(d_1=1),..., P(F+,F-|d_n=1)p(d_n=1)
# the 2nd until (n+1)nd element are the terms for d_1=1,...,d_n=1 and for all the elements, normalization is done with p(F+,F-).
# posterior is the most important output: posterior = P(d_1=1|F+,F-),... P(d_n=1|F+,F-)
