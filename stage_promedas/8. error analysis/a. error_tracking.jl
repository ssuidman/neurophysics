include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl") 
include("../0. useful/run_one_time.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")


function repeat1(;n=10)
    """
    This function shows how eps(Float64) gives the smallest steps that can be taken (~1e-16)
    """
    y = 1.0
    z = 1.0
    println(z-1,"\t\t\t",y-1,"\t\t\t",eps(Float64))
    for i=1:n
        y = prevfloat(y)
        z = nextfloat(z)
        println(z-1,"\t",y-1,"\t",eps(Float64))
    end
end 
repeat1()


patient_cases_raw, data_alisa, previn1, pfmin1, pfminneg1, sens1, sensneg1, prev1, prevminneg1 = prepare_patient_data("case 1");
previn2, pfmin2, pfminneg2, sens2, sensneg2, prev2, prevminneg2 = quickscore_preparation(size(pfmin1,1),n_disease=size(pfmin1,2));
# pfplus_prod, P_joint_prod, posterior_prod, dt_prod = quickscore(previn1, pfmin1, pfminneg1,"prod");  
term1 = prod(1e-50 .+ (prod(pfmin1,dims=1) .* prevminneg1 .+ (1 .- prev1)),dims=2);
term1_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin1),dims=1) .* prevminneg1 .+ (1 .- prev1)),dims=2);
maximum(abs,term1 - term1_BF)

term2 = prod(1e-50 .+ (prod(pfmin2,dims=1) .* prevminneg2 .+ (1 .- prev2)),dims=2);
term2_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin2),dims=1) .* prevminneg2 .+ (1 .- prev2)),dims=2);
maximum(abs,term2 - term2_BF)

pfmin3, prevminneg3, prev3 = rand(size(pfmin1,1),size(pfmin1,2)), rand(size(prevminneg1,1),size(prevminneg1,2)), rand(size(prev1,1),size(prev1,2));
term3 = prod(1e-50 .+ (prod(pfmin3,dims=1) .* prevminneg3 .+ (1 .- prev3)),dims=2); 
term3_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin3),dims=1) .* prevminneg3 .+ (1 .- prev3)),dims=2); 
maximum(abs,term3 - term3_BF)
k = 10000;
pfmin4, prevminneg4, prev4 = 1 .- rand(size(pfmin1,1),size(pfmin1,2))/k, 1 .- rand(size(prevminneg1,1),size(prevminneg1,2))/k, 1 .- rand(size(prev1,1),size(prev1,2))/k;
term4 = prod(1e-50 .+ (prod(pfmin4,dims=1) .* prevminneg4 .+ (1 .- prev4)),dims=2); 
term4_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin4),dims=1) .* prevminneg4 .+ (1 .- prev4)),dims=2); 
maximum(abs,term4 - term4_BF)


# previn, pfmin, prevminneg, prev, pfplus, pfplus_BF = fill(1. - eps(Float64),size(previn)), fill(1. - eps(Float64),size(pfmin)), fill(1. - eps(Float64),size(prevminneg)), fill(1. - eps(Float64),size(prev)), zeros(n+1,1), zeros(n+1,1); 
# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(7,n_disease=773);
patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
denominators = 2.220603782519935525508225413674727089837893977516790520311692475798540222583655e-08
m,n = size(pfmin); 
pfplus, pfplus_BF = zeros(n+1,1), zeros(n+1,1); 
for i in ProgressBar(0:(2^m-1)) 
    v = digits(i,base=2,pad=m) 
    myset = findall(v.==1) 
    # myset = 1:m
    pfplus = pfplus .+ ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
    pfplus_BF = pfplus_BF .+ ((-1)^length(myset)) .* prod(1e-50 .+ (prod(BigFloat.(pfmin[myset, :]),dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
    
    posterior = pfplus[2:end,1] .* previn ./ denominators
    posterior_BF = pfplus_BF[2:end,1] .* previn ./ denominators
    max = maximum(abs,posterior) 
    if max > 1
        println(myset,"\t",Float64(max),"\t",Float64.(pfplus)[1:5])
        println(myset,"\t",Float64(max),"\t",Float64.(pfplus)[1:5])
    end
end

posterior = pfplus[2:end,1] .* previn ./ pfplus[1,1]
posterior_BF = pfplus_BF[2:end,1] .* previn ./ pfplus_BF[1,1]

posterior[369-2:369+2]
posterior_BF[369-2:369+2]

pfplus_BF[2:end,1] .* previn ./ denominators

# P_joint_BF = pfplus_BF[2:end,1] .* previn; 
# posterior_BF = P_joint_BF / pfplus_BF[1,1]; 
# P_joint = pfplus[2:end,1] .* previn; 
# posterior = P_joint / pfplus[1,1]; 

myset = [1,2,3];
(((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2))[2:end,1] .* previn 
# posterior - posterior_BF

t0 = prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
t1 = prod(pfmin,dims=1) .* (1 .+ prevminneg .* sum(1 ./ pfmin,dims=1)) .+ 1
t2 = prod(pfmin,dims=1) .* prevminneg .+ (1 .- prev)
t3 = t0.*sum(t1./t2,dims=2) 
t4 = t3 * 2^m 



for (i,post) in enumerate(Float64.(posterior))
    if Float64.(abs.(post.-posterior_BF[i])) > 0.0001
    # if post > 1 || post < 0
        println("############### error ###############")
        println(i,"\t",round.(post,sigdigits=1),"\t\t",Float64.(abs.(post.-posterior_BF[i])),"\t\t",pfplus[i],"\t\t",Float64.(pfplus_BF[i]))
        # println(i,"\t",round.(post,sigdigits=1),"\t\t",Float64.(abs.(post.-posterior_BF[i])),"\t\t",P_joint[i],"\t\t\t",Float64.(P_joint_BF[i]))
    else
        # println(i,"\t",round.(post,sigdigits=1),"\t\t",Float64.(abs.(post.-posterior_BF[i])),"\t\t",pfplus[i],"\t\t",Float64.(pfplus_BF[i]))
    end
end





include("../0. useful/run_one_time.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")
# previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");
# patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(20,n_disease=1000);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"trick BF thread"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128"); 

posterior - posterior_

m = 18
previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(m,n_disease=1000);
v_ = digits.(0:2^m-1,base=2,pad=m)
myset_ = Vector{Vector{Int64}}(undef, 2^m)
pfmin_prod = Array{Float64,3}(undef,(2^m,1,size(pfmin,2)))
for i in ProgressBar(1:2^m)
    myset_[i] = findall(v_[i].==1)
    pfmin_prod[i,:,:] = prod(pfmin[myset_[i],:],dims=1)
end
myset = myset_[100]
println("Sizes:"," myset: ",size(myset)," pfmin: ",size(pfmin)," prevminneg: ",size(prevminneg)," prev: ",size(prev))



println("a1");  @btime a1 = pfmin[myset, :]; a1 = pfmin[myset, :]; 
println("a2");  @btime a2 = prod(a1, dims=1); a2 = prod(a1, dims=1);
println("b");   @btime b = prevminneg; b = prevminneg; 
println("1 - prev");   @btime c = 1 .- prev; c = 1 .- prev; # This step takes long (161 μs)
println("prod(pfmin) * prevminneg");   @btime d = a2 .* b; d = a2 .* b; # This step takes long (162 μs)
println("prof(pfmin)*prevminneg + (1-prev)");   @btime e = d .+ c;    e = d .+ c; # This step takes long (258 μs)
println("prod(inner_term)");   @btime f = prod(e,dims=2); f = prod(e,dims=2); # This step takes long (145 μs)
println("g");   @btime g = ((-1)^length(myset)); g = ((-1)^length(myset)); 
println("h");   @btime h = g .* f; h = g .* f; 




##### HIER VERDER!!! #####
include("../0. useful/run_one_time.jl"); 
include("../1. random m=1-22/a. quickscore_preparation.jl"); 
include("../0. useful/quickscore_algorithm.jl"); 
include("../5. C++ test/a. cpp_preparation_functions.jl"); 
# previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");
patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"prod BF thread"); println(typeof(posterior)); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick BF"); println(typeof(posterior)); 

prev1 = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); 
prev2 = BigFloat.(previn)' .+ [zeros(n) I(n)]' .* (1 .- BigFloat.(previn)'); 



# Run one time 
m,n = size(pfmin); 
myset = [1,2,5,7]; 
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); 
previn, pfminneg, pfmin, prev = [BigFloat.(x) for x in [previn, pfminneg, pfmin, prev]]; 

prevminneg = prev .* pfminneg'; 
prod_pfmin = prod(pfmin[myset,:],dims=1); 
x_old = prod(pfmin[myset,:],dims=1) .* prevminneg .+ (1 .- prev); 
x_new = prod_pfmin .* previn' .* pfminneg' .+ (1 .- previn') .+ (prod_pfmin .* pfminneg' .- 1) .* (1 .- previn') .* [zeros(n) I(n)]'; 

previn_pfminneg = previn' .* pfminneg'; 
one_min_previn = 1 .- previn'; 
pfminneg_min_previn = pfminneg' .- previn'; 

y_old = ((-1)^length(myset)) .* prod(prod(pfmin[myset,:],dims=1) .* prevminneg .+ (1 .- prev),dims=2); 
A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
B = (prod_pfmin .* pfminneg' .- 1) .* (1 .- previn'); 
y_new = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']; 
sum(abs, x_new-x_old) 
sum(abs, y_new-y_old)




# Run quickscore algorithm 
m,n = size(pfmin); # done
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn'); # done
previn, pfminneg, pfmin, prev = [BigFloat.(x) for x in [previn, pfminneg, pfmin, prev]]; 
prevminneg = prev .* pfminneg'; # done

# previn, pfminneg, pfmin = [BigFloat.(x) for x in [previn, pfminneg, pfmin]]; # This line makes it different 

pfplus_matrix_new = zeros(BigFloat,2^m,n+1,1); 
pfplus_matrix_old = zeros(BigFloat,2^m,n+1,1); 

previn_pfminneg = previn' .* pfminneg'; 
one_min_previn = 1 .- previn'; 
pfminneg_min_previn = pfminneg' .- previn'; 

# s = 0
# for i in ProgressBar(0:2^m-1)
for i in 0:2^m-1
    v = digits(i,base=2,pad=m)
    myset = findall(v .== 1)
    prod_pfmin = prod(pfmin[myset,:],dims=1) 
    A = prod_pfmin .* previn_pfminneg .+ one_min_previn
    B = (prod_pfmin .* pfminneg' .- 1) .* (1 .- previn')
    pfplus_matrix_new[i+1,:,:] = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']
    pfplus_matrix_old[i+1,:,:] = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2)
    # result_old = ((-1)^length(myset)) .* prod(prod(pfmin[myset,:],dims=1) .* prevminneg .+ prev_min_1,dims=2); 
    println(maximum(abs,pfplus_matrix_new[i+1,:,:] .- pfplus_matrix_old[i+1,:,:]))
    # s = s .+ maximum(abs,result_new .- result_old)
end
sum(abs,pfplus_matrix_new .- pfplus_matrix_old)




# Directly from the quickscore function 
m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
prev = previn' .+ [zeros(n) I(n)]' .* (1 .- previn')
prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)

pfplus_matrix_new = zeros(Float64,2^m,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 
pfplus_matrix_old = zeros(Float64,2^m,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 

previn, pfmin, prevminneg, prev, pfplus, pfplus_matrix_new, pfplus_matrix_old  = [BigFloat.(x) for x in [previn, pfmin, prevminneg, prev, pfplus, pfplus_matrix_new, pfplus_matrix_old]]; 

previn_pfminneg = previn' .* pfminneg' 
one_min_previn = 1 .- previn'
pfminneg_min_previn = pfminneg' .- previn'

dt = @elapsed begin
    for i in ProgressBar(0:(2^m-1)) # iterate over 2^m possibilities 
        v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
        myset = findall(v.==1) # find places of the 1-elements
        pfplus_matrix_old[i+1,:,:] = ((-1)^length(myset)) .* prod(prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev),dims=2)
        prod_pfmin = prod(pfmin[myset,:],dims=1) 
        A = prod_pfmin .* previn_pfminneg .+ one_min_previn
        B = (prod_pfmin .* pfminneg' .- 1) .* (1 .- previn')
        pfplus_matrix_new[i+1,:,:] = ((-1)^length(myset)) .* prod(A) .* [1 ; (A .+ B)' ./ A']
    end
end
pfplus_matrix_old .- pfplus_matrix_new





