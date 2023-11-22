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
previn2, pfmin2, pfminneg2, sens2, sensneg2, prevminneg2 = quickscore_preparation(size(pfmin1,1),n_disease=size(pfmin1,2));
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

pfmin3, prevminneg3, prev3 = rand(size(pfmin1,1),size(pfmin1,2)), rand(size(prevminneg1,1),size(prevminneg1,2)), rand(size(prev1,1),size(prev1,2));
term3 = prod(1e-50 .+ (prod(pfmin3,dims=1) .* prevminneg3 .+ (1 .- prev3)),dims=2); 
term3_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin3),dims=1) .* prevminneg3 .+ (1 .- prev3)),dims=2); 
maximum(abs,term3 - term3_BF)

k = 10000;
# k = 1000;
pfmin4, prevminneg4, prev4 = 1 .- rand(size(pfmin1,1),size(pfmin1,2))/k, 1 .- rand(size(prevminneg1,1),size(prevminneg1,2))/k, 1 .- rand(size(prev1,1),size(prev1,2))/k;
term4 = prod(1e-50 .+ (prod(pfmin4,dims=1) .* prevminneg4 .+ (1 .- prev4)),dims=2); 
term4_BF = prod(1e-50 .+ (prod(BigFloat.(pfmin4),dims=1) .* prevminneg4 .+ (1 .- prev4)),dims=2); 
Float64.(maximum(abs,term4 - term4_BF))

1 .- rand()/k

