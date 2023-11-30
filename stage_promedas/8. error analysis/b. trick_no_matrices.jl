include("../0. useful/packages.jl") 
include("../0. useful/run_one_time.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(9,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick BF"); 
# pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread"); 
# posterior - posterior_



##############################################################################################################################
############## Decompose all the steps of the quickscore algorithm, so you can look what takes the longest time ##############
##############################################################################################################################
println("pfmin[myset,:]");  @btime a1 = pfmin[myset, :]; a1 = pfmin[myset, :]; 
println("prod( pfmin[myset,:] )");  @btime a2 = prod(a1, dims=1); a2 = prod(a1, dims=1);
println("1 - prev");   @btime c = 1 .- prev; c = 1 .- prev; # This step takes long (161 μs)
println("prod(pfmin) * prevminneg"); b = prevminneg; @btime d = a2 .* b; d = a2 .* b; # This step takes long (162 μs)
println("prof(pfmin)*prevminneg + (1-prev)");   @btime e = d .+ c;    e = d .+ c; # This step takes long (258 μs)
println("prod(inner_term)");   @btime f = prod(e,dims=2); f = prod(e,dims=2); # This step takes long (145 μs)
println("(-1)^length(myset)");   @btime g = ((-1)^length(myset)); g = ((-1)^length(myset)); 
println("(-1)^length(myset) * prod(inner_term)");   @btime h = g .* f; h = g .* f; 
println("Conclusion:\n\tTerms that take the longest time are \n\t\t'1 - prev'\n\t\t'prod(pfmin) * prevminneg'\n\t\t'prod(pfmin)*prevminneg + (1-prev)'\n\t\t'prod(inner_term)'\n\tThe reason is that they are (n+1,n) matrices.")





###############################################################################################################################
########################################## Let's see if we can decompose the matrix  ##########################################
###############################################################################################################################

########################### We start with variables in BigFloat: ###########################
previn, pfminneg, pfmin = [BigFloat.(x) for x in [previn, pfminneg, pfmin]]; 
########################### Then we calculate prev: ###########################
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); 
########################### Then for prevminneg we fill in prev: ###########################
prevminneg  = prev .* pfminneg'; 
#           = (previn' .+ [zeros(n) I(n)]' .* (1 .- previn')) .* pfminneg' 
#           = previn' .* pfminneg'   .+     (1 .- previn) .* pfminneg' .* [zeros(n) I(n)]'
########################### Then for the inner_term we fill in prevminneg, prev and recollect terms: ###########################
inner_term  = prod(pfmin[myset,:],dims=1) .* prevminneg .+ (1 .- prev); 
#           = prod(pfmin[myset,:],dims=1) .* (previn' .* pfminneg'  .+  (1 .- previn) .* pfminneg' .* [zeros(n) I(n)]')     .+      (1   .-   (previn' .+ [zeros(n) I(n)]' .* (1 .- previn')) )
#           = prod(pfmin[myset,:],dims=1) .* previn' .* pfminneg' .+ (1 .- previn')     .+      (prod(pfmin[myset,:],dims=1) .* pfminneg' .- 1) .* (1 .- previn') .* [zeros(n) I(n)]'
########################### Precalculate some variables before running quickscore   ###########################
previn_pfminneg = previn' .* pfminneg'; 
one_min_previn = 1 .- previn'; 
pfminneg_min_previn = pfminneg' .- previn'; 
########################### Rewrite 'inter_term': ###########################
#           = prod(pfmin[myset,:],dims=1) .* previn_pfminneg .+ one_min_previn  .+  (prod(pfmin[myset,:],dims=1) .* pfminneg' .- 1) .* one_min_previn .* [zeros(n) I(n)]'
########################### Calculate values during quickscore: ###########################
i = 13; 
v = digits(i,base=2,pad=m); 
myset = findall(v .== 1);
prod_pfmin = prod(pfmin[myset,:],dims=1); 
A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
B = (prod_pfmin .* pfminneg' .- 1) .* (1 .- previn'); 
########################### Such that 'inter_term' consists of a vector and a diagonal: ###########################
#           = A .+ B .* [zeros(n) I(n)]'
println("A:")
pretty_table(Float64.(A))
println("A+B:")
pretty_table(Float64.(A .+ B))
println("A + B * diagn")
pretty_table(Float64.(A .+ B .* [zeros(n) I(n)]'))
########################### Calculate first-row product, then only flip one term for each row, multiply by +/-1 ###########################
term_new = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; 
term_old = ((-1)^length(myset)) .* prod(prod(pfmin[myset,:],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
########################### Compare the new term with the old term ###########################
println("Maximum difference in 1 term to be summed into pfplus between old and new method is\n\t",maximum(abs,term_new .- term_old))
println("Conclusion:\n\tThe (n+1,n) matrices can be decomposed in 2 (1,n) matrices of which the product can be taken efficiently")




