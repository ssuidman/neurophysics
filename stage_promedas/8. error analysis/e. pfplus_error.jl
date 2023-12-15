include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")



patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)];
previn, pfmin, pfminneg = rand(n)/100, 1 .-rand(m,n)/100, rand(n)
pfplus, x, z = zeros(n+1,1), copy(previn)', copy(pfminneg)';
for (i,myset) in enumerate(myset_matrix)
    Y = prod(pfmin[myset,:],dims=1)
    # term = (-1)^length(myset).*prod(Y.*x.*z.+(1 .- x),dims=2).*[1 (Y.*z)./(Y.*z.*x.+(1 .- x))]
    # term = (-1)^length(myset).*prod(Y.*x.*z.+(1 .- x),dims=2)
    d = 1 ./1000
    term = (-1)^length(myset).*prod(1-d .+ Y.*d,dims=2)
    pfplus = pfplus .+ term'
    signs = (-1)^length(myset)==1 ? "+" : "-"
    # println(i-1,"\t",digits(i-1,base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:4],sigdigits=2)',"\t",round.(pfplus[1:4],sigdigits=2)')
    println(i-1,"\t",digits(i-1,base=2,pad=m),"\t",signs," ",round.(pfplus[1:5],sigdigits=3)')
end

prod_pfmin = prod(pfmin[[],:],dims=1)
term_0 = (-1)^length(myset).*prod(prod_pfmin.*previn'.*pfminneg'.+(1 .-previn'),dims=2).*[1 (prod_pfmin.*pfminneg')./(prod_pfmin.*pfminneg'.*previn'.+(1 .-previn'))]
prod_pfmin = prod(pfmin[[1:m...],:],dims=1)
term_last = (-1)^length(myset).*prod(prod_pfmin.*previn'.*pfminneg'.+(1 .-previn'),dims=2).*[1 (prod_pfmin.*pfminneg')./(prod_pfmin.*pfminneg'.*previn'.+(1 .-previn'))]
