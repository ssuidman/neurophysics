include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(20,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick Fl128"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 dict"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread"); 
# posterior .- posterior_



################################################################################################################
################################ Look at the creation time of prod_pfmin_matrix ################################
################################################################################################################
previn, pfmin, pfminneg = [Float128.(x) for x in [previn, pfmin, pfminneg]];  
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); 
prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n); # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
previn_pfminneg = previn' .* pfminneg'; 
one_min_previn = 1 .- previn'; 
pfminneg_min_previn = pfminneg' .- previn'; 
pfplus_matrix = zeros(typeof(previn[1]),2^m,n+1,1); 
prod_pfmin_matrix = zeros(typeof(previn[1]),2^m,1,n); 
v_matrix = permutedims(hcat(digits.(0:2^m-1,base=2,pad=m)...)); 

dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements 
        prod_pfmin_matrix[i+1,:,:] = prod(pfmin[myset,:],dims=1); 
    end 
end; 
println("Running time for the creation of prod_pfmin_matrix:\n\t",dt," sec")

dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements 
        prod_pfmin = prod(pfmin[myset,:],dims=1); 
        A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
        B = (prod_pfmin .* pfminneg' .- 1) .* one_min_previn; 
        pfplus_matrix[i+1,:,:] = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; 
    end 
end; 
println("Running time quickscore $(typeof(previn[1])) normal:\n\t",dt," sec")

dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements 
        prod_pfmin = prod_pfmin_matrix[i+1,:,:]
        A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
        B = (prod_pfmin .* pfminneg' .- 1) .* one_min_previn; 
        pfplus_matrix[i+1,:,:] = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; 
    end 
end; 
println("Running time quickscore $(typeof(previn[1])) prod_pfmin_matrix:\n\t",dt," sec")
println("Conclusion: For Float128, creating prod_pfmin is taking half the time. If we can get O(m^2) instead of O(2^m), this is beneficial for higher m.")
println("Therefore we should create a faster for-loop for prod_pfmin_matrix")



################################################################################################################
################## The problem is now to speed up the for-loop for creating prod_pfmin_matrix ##################
################################################################################################################
float_type = typeof(previn[1]); 
prod_pfmin_matrix = zeros(float_type,2^m,1,n); 
dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements 
        prod_pfmin_matrix[i+1,:,:] = prod(pfmin[myset,:],dims=1); 
    end 
    println("Running time $float_type normal:")
end

myset_matrix = collect(powerset([1:m...])); 
# prod_pfmin_dict = Dict(myset_matrix[i] => prod(pfmin[myset_matrix[i],:],dims=1) for i in 1:m+1); 
prod_pfmin_dict = Dict{Vector{Int64},Matrix{typeof(previn[1])}}(myset_matrix[i] => prod(pfmin[myset_matrix[i],:],dims=1) for i in 1:m+1); 
@elapsed dt = begin
    for myset in ProgressBar(myset_matrix[m+2:end])
        prod_pfmin = prod_pfmin_dict[myset[1:end-1]] .* prod_pfmin_dict[myset[[end]]]  
        push!(prod_pfmin_dict,myset => prod_pfmin)
    end
    println("Running time $float_type dict:")
end

myset = myset_matrix[100]
prod_pfmin_dict[myset]
prod(pfmin[myset,:],dims=1)

myset_matrix = collect(powerset([1:21...]))



