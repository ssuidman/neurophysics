# Wim Wiegerinck version 8/9/2023 -- adapted Sam Suidman 25/9/2023 
# inputs:
# previn               prevalences
# pfmin = 1 - sens     (for positive tests): these are Heckerman's p(f-|d)
# pfminneg             (products of 1-sens for negatives, to be absorbed by prevs): these are hekkermans Prod_{F-} (p(f-|d))
# outputs:
# posterior =  [P(d_1=1|F+,F-),... P(d_n=1|F+,F-)]
# pfpluspd = [p(F+,F-), P(F+,F-|d_1=1)p(d_1=1),..., P(F+,F-|d_n=1)p(d_n=1)]
# pfplus = [ p(F+,F-), P(F+,F-|d_1=1),..., P(F+,F-|d_n=1) ]
using ProgressBars

function quickscore(previn, pfmin, pfminneg,binary)
    m, n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevend = vcat(1,previn) # prevend are the prevalences 1, P(d_1=1),... P(d_n=1) and is needed to multiply at the end
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    
    pfplus = zeros(n+1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
    myset = 0
    t1 = time()
    if binary
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
    else
        ready = false 
        v = zeros(Int,m) 
        while !ready
            myset = findall(v.==1) # find places of the 1-elements
            if length(myset)==0
                pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ prevminneg .+ (1 .- prev)), dims=2))
            elseif length(myset)==1
                pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev)), dims=2))
            else
                pfplus .+= ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
            end
    
            ready = true 
            for k = 1:m
                v[k] += 1
                if v[k] <= 1
                    ready = false
                    break
                end
                if k > 1
                    v[k] = v[k - 1]
                else
                    v[k] = 0
                end
            end
        end        
    end     
    t2 = time()
    dt = t2 - t1 
    
    pfpluspd = pfplus.*prevend 
    posterior = pfpluspd[2:end] / pfpluspd[1] 
    
    return pfplus, posterior, dt 
end 

# final results
# pfplus is p(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1),
# pfpluspd is numerator in heckerman (12)  p(F+,F-),P(F+,F-|d_1=1)p(d_1=1),..., P(F+,F-|d_n=1)p(d_n=1)
# the 2nd until (n+1)nd element are the terms for d_1=1,...,d_n=1 and for all the elements, normalization is done with p(F+,F-).
# posterior is the most important output: posterior = P(d_1=1|F+,F-),... P(d_n=1|F+,F-)
