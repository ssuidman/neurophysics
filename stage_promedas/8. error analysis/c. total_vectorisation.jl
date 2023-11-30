include("../0. useful/packages.jl") 
include("../0. useful/run_one_time.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 3"); myset = [1,2,3]; m,n=size(pfmin);
previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(18,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick BF "); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"prod BF"); 
maximum(abs,posterior - posterior_)


prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn')
prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
previn_pfminneg = previn' .* pfminneg' 
one_min_previn = 1 .- previn'
pfminneg_min_previn = pfminneg' .- previn'

@btime prod_pfmin = prod(pfmin[myset,:],dims=1); 
@btime A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
@btime B = (prod_pfmin .* pfminneg' .- 1) .* one_min_previn;  
@btime term = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A'];



prod_pfmin_matrix = Array{Float64,3}(undef,(2^m,1,size(pfmin,2))); 
v_matrix = digits.(0:2^m-1,base=2,pad=m); 
myset_matrix = Vector{Vector{Int64}}(undef, 2^m); 
dt = @elapsed begin
    for i in ProgressBar(0:2^m-1)
        myset_matrix[i+1] = findall(v_matrix[i+1].==1)
        prod_pfmin_matrix[i+1,:,:] = prod(pfmin[myset_matrix[i+1],:],dims=1)
    end
end; 
println("Running time: ",dt)

previn_pfminneg = previn' .* pfminneg'; 
one_min_previn = 1 .- previn'; 
pfminneg_min_previn = pfminneg' .- previn'; 
pfplus_matrix = zeros(Float64,2^m,n+1,1) ; 
dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        # myset = findall(v.==1) # find places of the 1-elements 
        # myset = findall(v[i].==1) # find places of the 1-elements 
        prod_pfmin = prod_pfmin_matrix[i+1,:,:]; 
        # prod_pfmin = prod(pfmin[myset,:],dims=1); 
        A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
        B = (prod_pfmin .* pfminneg' .- 1) .* one_min_previn; 
        pfplus_matrix[i+1,:,:] = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; 
    end 
end; 
println("Running time: ",dt)





previn, pfminneg, pfmin = [BigFloat.(x) for x in [previn, pfminneg, pfmin]]; 
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); 
prevminneg  = prev .* pfminneg'; 

v_matrix = permutedims(hcat(digits.(0:2^m-1,base=2,pad=m)...))
1 .- reshape(v_matrix,size(v_matrix)...,1)
reshape(v_matrix,size(v_matrix)...,1)

reshape(v_matrix,size(v_matrix)...,1) .* reshape(pfmin,1,size(pfmin)...)

reshape(pfmin,1,size(pfmin)...) 
prod_pfmin_matrix = prod(reshape(v_matrix,size(v_matrix)...,1) .* reshape(pfmin,1,size(pfmin)...) .+ (1 .- v_matrix), dims=2)[:,1,:]
A = prod_pfmin_matrix .* pfminneg' .* previn' .+ (1 .- previn')
B = (prod_pfmin_matrix .* pfminneg' .- 1) .* (1 .- previn')
pfplus_matrix = (-1).^sum(v_matrix,dims=2) .* prod(A,dims=2) .* [1 ; (A .+ B)' ./ A']'




for i=1:512
    idx = findall(digits(i-1,base=2,pad=m) .== 1)
    maxx = maximum(abs,x[i,:] .- prod(prod(pfmin[idx,:],dims=1) .* prevminneg .+ (1 .- prev),dims=2))
    println(i," ",maxx)
end
# This should be pfplus: 
maximum(abs,sum(x .* (-1).^sum(v_matrix,dims=2),dims=1)' .- pfplus)


pfplus_matrix

for i=1:512

    prod_pfmin_matrix[i]


A[1,:]' .+ B[1,:]' .* [zeros(n) I(n)]'
prod(A,dims=2)
A[1,:]

# A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
# B = (prod_pfmin .* pfminneg' .- 1) .* one_min_previn;
# term = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; # --> prod(A + B * I0), with 'I0' the diagonal. 



prod(prod_pfmin_matrix .* pfminneg',dims=2)


test[test .== 0.0] .= 1.0
prod(test,dims=2)[23+1,1,:]'
prod(pfmin[[1,2,3,5],:],dims=1)
digits(23,base=2,pad=10)

prod(pfmin[1,:],dims=2)
test[1,:]


digits(2^(3-1)+2^(8-1)+2^(6-1),base=2,pad=10)'








