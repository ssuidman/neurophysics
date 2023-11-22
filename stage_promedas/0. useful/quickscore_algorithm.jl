include("packages.jl")
include("useful_functions.jl")

function loop_term(myset,pfmin,prevminneg,prev,method)
    if      method in ["exp-sum-log"];                  term = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
    elseif  method in ["exp-sum-log prod1"];            term = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
    elseif  method in ["exp-sum-log Fl32"];             term = ((-1)^length(myset)) .* exp.(sum(log.(Float32(1e-50) .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
    elseif  method in ["exp-sum-log QM Fl128"];         term = ((-1)^length(myset)) .* exp.(sum(log.(Float128(1e-50) .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
    elseif  method in ["exp-sum-log BF(entire term)"];  term = BigFloat.(((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2)));
    elseif  method in ["prod"];                         term = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
    elseif  method in ["prod prod1"];                   term = ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2);  
    elseif  method in ["prod Fl32"];                    term = ((-1)^length(myset)) .* prod(Float32(1e-50) .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
    elseif  method in ["prod QM Fl128"];                   term = ((-1)^length(myset)) .* prod(Float128(1e-50) .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
    elseif  method in ["prod BF(entire term)"];         term = BigFloat.(((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2));
    elseif  method in ["prod BF","prod BF Fl64","prod BF Fl128"];   term = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(BigFloat.(pfmin)[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
    elseif  method in ["exp-sum-log BF"];               term = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(BigFloat.(pfmin)[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2)); 
    end
    return term
end

function quickscore(previn::Vector{Float64}, pfmin::Matrix{Float64}, pfminneg::Vector{Float64}, method::String)
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $m (postive tests)")
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    pfplus = zeros(Float64,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    
    if      occursin("Fl32",method);    previn, pfmin, prevminneg, prev, pfplus = Float32.(previn), Float32.(pfmin), Float32.(prevminneg), Float32.(prev), Float32.(pfplus)
    elseif  occursin("QM Fl128",method);previn, pfmin, prevminneg, prev, pfplus = Float128.(previn), Float128.(pfmin), Float128.(prevminneg), Float128.(prev), Float128.(pfplus)
    elseif  occursin("Fl64",method);    setprecision(BigFloat,53)
    elseif  occursin("BF Fl128",method);setprecision(BigFloat,113)
    else;                               setprecision(BigFloat,256)
    end
    
    if method!="MATLAB"
        println("Method: '$method'")
        dt = @elapsed begin
            for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
                v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
                myset = findall(v.==1) # find places of the 1-elements
                pfplus = pfplus .+ loop_term(myset,pfmin,prevminneg,prev,method)

                # term = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(BigFloat.(pfmin[myset, :]),dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 

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

    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("\nRunning time: $dt")
    return pfplus, P_joint, posterior, dt 
end 
# pfplus_, P_joint_, posterior_,dt_ = quickscore(previn, pfmin, pfminneg,"prod");  


# final results
# pfplus is p(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1),
# pfpluspd is numerator in heckerman (12)  p(F+,F-),P(F+,F-|d_1=1)p(d_1=1),..., P(F+,F-|d_n=1)p(d_n=1)
# the 2nd until (n+1)nd element are the terms for d_1=1,...,d_n=1 and for all the elements, normalization is done with p(F+,F-).
# posterior is the most important output: posterior = P(d_1=1|F+,F-),... P(d_n=1|F+,F-)
