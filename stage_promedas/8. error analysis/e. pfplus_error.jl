include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

# First look at an example of the actual case 
patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
myset_matrix = [powerset(1:m)...]
pfplus = zeros(n+1,1)
for (i,myset) in enumerate(myset_matrix)
    prod_pfmin = prod(pfmin[myset,:],dims=1)
    term = ((-1)^length(myset)) .* prod(prod_pfmin .* previn'.*pfminneg' .+ (1 .- previn')) .* [1 ; (prod_pfmin .* pfminneg')' ./ (prod_pfmin .* previn'.*pfminneg' .+ (1 .- previn'))']; # --> prod(A + (B-A) * I0), with I0=[zeros(n) I(n)]' the diagonal. 
    pfplus = pfplus .+ term'
    signs = (-1)^length(myset)==1 ? "+" : "-"
    println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(Float64.(pfplus)[1:5],sigdigits=3)')
    # println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:4],sigdigits=3)',"\t",round.(pfplus[1:4],sigdigits=3)')
end



# The problem can be brought back to this simple problem 
myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)];
d1, d2 = 1 ./1000000000000, 1 ./100 # changing d2 with one zero (factor 10) gives an extra 1e-7 (1e-m) in the end, the difference for changing d1 becomes smaller and smaller 
pfplus = 0 # zeroth element 
pfmin = BigFloat.(1 .- rand(m,n).*d2)
for (i,myset) in enumerate(myset_matrix)
    prod_pfmin = prod(pfmin[myset,:],dims=1)
    # term = (-1)^length(myset).*prod(1-d1 .+ prod_pfmin.*d1,dims=2)
    term = (-1)^length(myset).*(1-d1)^n.*(1 .+ d1/(1-d1).*(1+sum(prod(pfmin[[1,3,4],:],dims=1))))
    pfplus = pfplus .+ term[1]
    signs = (-1)^length(myset)==1 ? "+" : "-"
    println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(Float64.(pfplus),sigdigits=3)')
    # println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:4],sigdigits=3)',"\t",round.(pfplus[1:4],sigdigits=3)')
end
(1-d1)^(n-1)*d1
sum(prod(1 .-pfmin[[1:7...],:],dims=1))


# Now look at the term x_i1
prod(1-d1 .+ prod(pfmin[[1,3,4],:],dims=1).*d1,dims=2) # original
(1-d1)^n.*prod(1 .+ d1/(1-d1) .* prod(pfmin[[1,3,4],:],dims=1),dims=2) # written differently 
(1-d1)^n.*    (1 .+ d1/(1-d1) .* ( 1 + sum(prod(pfmin[[1,3,4],:],dims=1))) ) # first approximation 



