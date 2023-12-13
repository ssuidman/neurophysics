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



patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 3");
println("previn:\n\t $(round(mean(previn),sigdigits=3)) +/- $(round(std(previn),sigdigits=3))\npfmin:\n\t $(round(mean(pfmin),sigdigits=3)) +/- $(round(std(pfmin),sigdigits=3))\npfminneg:\n\t $(round(mean(pfminneg),sigdigits=3)) +/- $(round(std(pfminneg),sigdigits=3))")
gr()
plot(sort(previn),label="previn",legend=:right,color="blue",dpi=1200)
plot!(sort(pfminneg),label="pfminneg",color="green")
# for i=1:9
#     plot!(sort(pfmin,dims=2)'[:,i],label="pfmin $i",color="red",dpi=1200)
# end
plot!(mean(sort(pfmin,dims=2)',dims=2),label="pfmin avg",color="black")
plot!(sort(prod(pfmin,dims=1),dims=2)',label="prod_pfmin",color="black")
xlabel!("diagnoses sorted")
ylabel!("probability")



patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 3");
diff_1 = (prod(pfmin,dims=1) .- prod(BigFloat.(pfmin),dims=1)) ./ prod(BigFloat.(pfmin),dims=1)
diff_2 = (prod(pfmin,dims=1).*prevminneg .- prod(BigFloat.(pfmin),dims=1).*prevminneg) ./ (prod(BigFloat.(pfmin),dims=1).*prevminneg)
diff_3 = (prod(pfmin,dims=1).*prevminneg.+(1 .-prev) .- prod(BigFloat.(pfmin),dims=1).*prevminneg.+(1 .-prev)) ./ ((prod(BigFloat.(pfmin),dims=1).*prevminneg.+(1 .-prev)))
diff_4 = (prod(prod(pfmin,dims=1).*prevminneg.+(1 .-prev),dims=2) .- prod(prod(BigFloat.(pfmin),dims=1).*prevminneg.+(1 .-prev),dims=2)) ./ (prod(prod(BigFloat.(pfmin),dims=1).*prevminneg.+(1 .-prev),dims=2))

m,n = size(pfmin)
A = (prod(pfmin,dims=1).*pfminneg'.*previn'.+(1 .-previn')).*(1 .- [zeros(n) I(n)]') 
B = prod(pfmin,dims=1).*pfminneg'.*[zeros(n) I(n)]'
A .+ B

x = rand(m,n)/100
y = 1 .- rand(m,n)/100
z = 1 .- rand(m,n)/100
prod(y,dims=1).*sum(1 ./y,dims=1).*eps(Float64)
abs.(prod(y,dims=1) .- prod(BigFloat.(y),dims=1))

a = rand(10).*100
a .- BigFloat.(a)
prod(a) .- prod(BigFloat.(a))
a


patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 1"); myset = [1,2,3]; m,n=size(pfmin);
# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(10,n_disease=1000); myset = [1,2,3]; m,n=size(pfmin);
pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore(previn, pfmin, pfminneg,"prod linear"); 
pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick BF linear"); 
# maximum(abs,posterior .- posterior_)
# hcat(posterior,posterior_)
myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)];

hcat([((-1)^length(myset_matrix[i])) .* prod(prod(pfmin[myset_matrix[i],:],dims=1) .* prevminneg .+ (1 .- prev),dims=2)[:,1] for i=1:2^m]...)'
prod(prod(pfmin[myset_matrix[4],:],dims=1) .* prevminneg .+ (1 .- prev),dims=2)

x = hcat([((-1)^length(myset_matrix[i])) .* prod(prod(pfmin[myset_matrix[i],:],dims=1) .* ones(size(prevminneg)),dims=2)[:,1] for i=1:2^m]...)'

sum(x,dims=1)
prod(pfmin[myset_matrix[5],:],dims=1)


pfplus_matrix_[:,:,1]

(pfplus_matrix[:,:,1] .- pfplus_matrix_[:,:,1]) ./ pfplus_matrix_[:,:,1]
pfplus_[:,1]'
((pfplus[:,1] .- pfplus_[:,1]) ./ pfplus_[:,1])'
# idx_zero = unique(hcat(collect.(Tuple.(findall(pfplus_matrix_[:,:,1].==0.0)))...)[2,:])[1]
# pfplus_matrix_[:,idx_zero-4:idx_zero+4,1]


myset = [1]
m,n = size(pfmin)
prod(prod(pfmin[myset,:],dims=1).*prevminneg.+(1 .-prev),dims=2)'
A = prod(pfmin[myset,:],dims=1) .* previn'.*pfminneg' .+ (1 .-previn'); 
B = prod(pfmin[myset,:],dims=1).*pfminneg'; 
# C = prod(A .+ (B .- A) .* [zeros(n) I(n)]',dims=2)'
C = prod(A) .* [1 ; B' ./ A']'
[1 ; B' ./ A']'


pfplus, x, z = zeros(n+1,1), copy(previn'), copy(pfminneg');
for (i,myset) in enumerate(myset_matrix)
    # prod_pfmin = prod(pfmin[myset,:],dims=1)
    Y = prod(pfmin[myset,:],dims=1)
    
    # term = (-1)^length(myset).*prod(prod_pfmin.*previn'.*pfminneg'.+(1 .-previn'),dims=2).*[1 (prod_pfmin.*pfminneg')./(prod_pfmin.*pfminneg'.*previn'.+(1 .-previn'))]
    # term = (-1)^length(myset).*[1 prod_pfmin]
    term = (-1)^length(myset).*prod(Y.*x.*z.+(1 .- x),dims=2).*[1 (Y.*z)./(Y.*z.*x.+(1 .- x))]
    
    pfplus = pfplus .+ term'
    signs = (-1)^length(myset)==1 ? "+" : "-"
    println(i-1,"\t",digits(i-1,base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:5],sigdigits=3)',"\t",round.(pfplus[1:5],sigdigits=3)')
    # println(i-1,"\t",digits(i-1,base=2,pad=m),"\t",round.(term[1:5]',sigdigits=3)')
    # println(i-1,"\t",round.(term'[1:5].-term_0'[1:5],sigdigits=3)')
    # println(i-1,"\t",round.(term'[1:5].-term_last'[1:5],sigdigits=3)')
end

prod_pfmin = prod(pfmin[[],:],dims=1)
term_0 = (-1)^length(myset).*prod(prod_pfmin.*previn'.*pfminneg'.+(1 .-previn'),dims=2).*[1 (prod_pfmin.*pfminneg')./(prod_pfmin.*pfminneg'.*previn'.+(1 .-previn'))]
prod_pfmin = prod(pfmin[[1:m...],:],dims=1)
term_last = (-1)^length(myset).*prod(prod_pfmin.*previn'.*pfminneg'.+(1 .-previn'),dims=2).*[1 (prod_pfmin.*pfminneg')./(prod_pfmin.*pfminneg'.*previn'.+(1 .-previn'))]


