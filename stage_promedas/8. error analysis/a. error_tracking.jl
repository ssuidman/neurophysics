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




myset = [1,2,3,4,5]
maximum(prod(pfmin[myset, :],dims=1))
maximum(prod(BigFloat.(pfmin[myset, :]),dims=1))




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





x1 = 1 ./ eps(Float64)^2
# x2 = BigFloat(1.1)^2e-16
# x3 = x2 - x1

significand(x1)
exponent(x1)    


# include("../0. useful/run_one_time.jl")
# previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");

######## HIER VERDER!!!!! #########
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")
patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1");
# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(9,n_disease=773);
pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"exp-sum-log BF",threading=true);  




Profile.init()
@profile pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"prod QM Fl128",threading=true);  
Profile.print()



