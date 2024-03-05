include("../0. useful/packages.jl") 
include("../0. useful/quickscore_algorithm.jl")
include("../1. random m=1-22/a. quickscore_preparation.jl")
include("../5. C++ test/a. cpp_preparation_functions.jl")

# First look at an example of the actual case 
patient_cases_raw, data_alisa, previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = prepare_patient_data("case 3"); myset = [1,2,3]; m,n=size(pfmin);
myset_matrix = [powerset(1:m)...]
myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)]
# previn = rand(n)./10
pfmin = rand(size(pfmin)...)./10
# pfminneg .=1
pfplus = zeros(n+1,1)
QQ = previn.*pfminneg #.+ rand(n)./1000
for (i,myset) in enumerate(myset_matrix)
    prod_pfmin = BigFloat.(prod(pfmin[myset,1:n],dims=1))
    # term = ((-1)^length(myset)) .* prod(prod_pfmin .* previn[1:n]'.*pfminneg[1:n]' .+ (1 .- previn[1:n]')) .* [1 ; (prod_pfmin .* pfminneg[1:n]')' ./ (prod_pfmin .* previn[1:n]'.*pfminneg[1:n]' .+ (1 .- previn[1:n]'))']; # --> prod(A + (B-A) * I0), with I0=[zeros(n) I(n)]' the diagonal. 
    # term = ((-1)^length(myset)) .* prod(1 .- previn[1:n]'.*(1 .- prod_pfmin .* pfminneg[1:n]')) .* [1 ; (prod_pfmin .* pfminneg[1:n]')' ./ (prod_pfmin .* previn[1:n]'.*pfminneg[1:n]' .+ (1 .- previn[1:n]'))']; # --> prod(A + (B-A) * I0), with I0=[zeros(n) I(n)]' the diagonal. 
    term = ((-1)^length(myset)) .* prod(1 .- previn[1:n]' .+ QQ'.*prod_pfmin ) .* [1 ; (prod_pfmin .* pfminneg[1:n]')' ./ (prod_pfmin .* previn[1:n]'.*pfminneg[1:n]' .+ (1 .- previn[1:n]'))']; # --> prod(A + (B-A) * I0), with I0=[zeros(n) I(n)]' the diagonal. 

    pfplus = pfplus .+ term
    signs = (-1)^length(myset)==1 ? "+" : "-"
    println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(Float64.(pfplus)[1:5],sigdigits=3)')
    # println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:4],sigdigits=3)',"\t",round.(pfplus[1:4],sigdigits=3)')
end

# rand(Beta(1,1000),10)
# plot(rand(fit(Beta,QQ),n))
Q1 =  prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[],1:10],dims=1) .* pfminneg[1:10]'))
Q2 = -prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[1],1:10],dims=1) .* pfminneg[1:10]'))
Q3 = -prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[2],1:10],dims=1) .* pfminneg[1:10]'))
Q4 =  prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[1,2],1:10],dims=1) .* pfminneg[1:10]'))
Q5 = -prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[3],1:10],dims=1) .* pfminneg[1:10]'))
Q6 =  prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[1,3],1:10],dims=1) .* pfminneg[1:10]'))
Q7 =  prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[2,3],1:10],dims=1) .* pfminneg[1:10]'))
Q8 = -prod(1 .- previn[1:10]'.*(1 .- prod(pfmin[[1,2,3],1:10],dims=1) .* pfminneg[1:10]'))
Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8
Q2 = prod(1 .- previn'.*(1 .- prod(pfmin[[1],:],dims=1) .* pfminneg'))
1 .- sum(previn'.*(1 .- pfmin[[1],:].*pfminneg')) + sum(prod(reduce(vcat,[previn[j]'.*(1 .- pfmin[[1],j].*pfminneg[j]') for j in combinations(1:n,2)]),dims=2))


pfmin[[1],:]




sum(previn'.*(1 .- pfminneg') .+ pfminneg'.*prod(1 .- pfmin[[5],:],dims=1))
1 .- pfmin
prod(pfmin[[1,2],:],dims=1)

# The problem can be brought back to this simple problem 
m,n = 2,2
myset_matrix = [findall(v.==1) for v in digits.(0:2^m-1,base=2,pad=m)];
d1, d2 = 1 ./1000, 1 ./100 # changing d2 with one zero (factor 10) gives an extra 1e-7 (1e-m) in the end, the difference for changing d1 becomes smaller and smaller 
pfplus = 0 # zeroth element 
# m,n = 7,773
pfmin = BigFloat.(1 .- rand(m,n).*d2)
for (i,myset) in enumerate(myset_matrix)
    prod_pfmin = prod(pfmin[myset,:],dims=1)
    term = (-1)^length(myset).*prod(1-d1 .+ prod_pfmin.*d1,dims=2)
    # term = (-1)^length(myset).*(1-d1)^n.*(1 .+ d1/(1-d1).*(1+sum(prod_pfmin)))
    pfplus = pfplus .+ term[1]
    signs = (-1)^length(myset)==1 ? "+" : "-"
    println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(Float64.(pfplus),sigdigits=3)')
    # println(i-1,"\t",digits(Int(sum(2 .^myset./2)),base=2,pad=m),"\t",signs," ",round.(abs.(term)[1:4],sigdigits=3)',"\t",round.(pfplus[1:4],sigdigits=3)')
end
pfplus
d1.*sum(prod(1 .-pfmin,dims=1))
# prod(pfmin[[2],:],dims=1).*pfminneg'.*previn'

# S0 = prod(1-d1 .+ prod(pfmin[[1,2,3],:],dims=1).*d1,dims=2)
S1a = prod(1 .-d1.*(1 .-prod(pfmin[[1],:],dims=1)),dims=2)
S1b = prod(1 .-d1.*(1 .-prod(pfmin[[2],:],dims=1)),dims=2)
S1c = prod(1 .-d1.*(1 .-prod(pfmin[[3],:],dims=1)),dims=2)
S2a = prod(1 .-d1.*(1 .-prod(pfmin[[1,2],:],dims=1)),dims=2)
S2b = prod(1 .-d1.*(1 .-prod(pfmin[[1,3],:],dims=1)),dims=2)
S2c = prod(1 .-d1.*(1 .-prod(pfmin[[2,3],:],dims=1)),dims=2)
S3  = prod(1 .-d1.*(1 .-prod(pfmin[[1,2,3],:],dims=1)),dims=2)
S = 1 .- S1a .- S1b .- S1c .+ S2a .+ S2b .+ S2c .- S3

prod(1 .-d1.*(1 .-prod(pfmin[[1],:],dims=1)),dims=2)[1,1]
1 .- d1.*(n.-sum(pfmin[[1],:])) .+ d2^2 .*(1/2*n*(n-1) .- (n-1).*sum(pfmin[[1],:]) .+ sum(prod(reduce(hcat,[combinations(pfmin[[1],:],2)...]),dims=1)))
1
d1.*(n.-sum(pfmin[[1],:]))
d2^2 .*(1/2*n*(n-1) .- (n-1).*sum(pfmin[[1],:]) .+ sum(prod(reduce(hcat,[combinations(pfmin[[1],:],2)...]),dims=1)))

x,a,b,c = 0,0,0,0
for (i,j) in ProgressBar(multiset_combinations([1:n...],2))
    x += (1-pfmin[1,i])*(1-pfmin[1,j])
    a += 1 # = binomial(n,2)
    b -= (pfmin[1,i] + pfmin[1,j]) # = -(n-1)*sum(pfmin[[1],:])
    c += pfmin[1,i]*pfmin[1,j] # = all combinations of 2 pfmin[[1],;]-terms --> sum(prod(reduce(hcat,[combinations(pfmin[[1],:],2)...]),dims=1))
end
y = binomial(n,2) -(n-1)*sum(pfmin[[1],:]) + sum(prod(reduce(hcat,[combinations(pfmin[[1],:],2)...]),dims=1))
x - y

expr_d1 = sum(prod(1 .-x,dims=2))
expr_d2 = (n-1)*sum(prod(1 .-x,dims=2)) - sum([prod(1 .-prod(x[i,1:end],dims=1)) for i in combinations(1:n,2)])

# expand(expr - expr_n_m)
mean(previn) # = d_1
std(previn) # = e_i
Float64(d1*subs(expr_d1,subs_dict))
Float64(d1^2*subs(expr_d2,subs_dict))
Float64((d1*subs(expr_d1,subs_dict) - d1^2*subs(expr_d2,subs_dict)))
Float64(d1*subs(expr_d1,subs_dict)/pfplus)
Float64((d1*subs(expr_d1,subs_dict) - d1^2*subs(expr_d2,subs_dict))/pfplus)

expr = 0; for i in ProgressBar(combinations(1:n,2))
for j in combinations(1:m)
    expr += (-1)^(length(j)+1)*prod(1 .- prod(x[i,j],dims=2))
end
end; 

# expr_x_line_1 = sum([prod(1 .- prod(x[i,[1]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]    # (1-x[1,1])*(1-x[2,1])+(1-x[1,1])*(1-x[3,1])+(1-x[2,1])*(1-x[3,1])
# expr_x_line_2 = sum([prod(1 .- prod(x[i,[2]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]    # (1-x[1,2])*(1-x[2,2])+(1-x[1,2])*(1-x[3,2])+(1-x[2,2])*(1-x[3,2])
# expr_x_line_3 = sum([prod(1 .- prod(x[i,[3]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]    # (1-x[1,3])*(1-x[2,3])+(1-x[1,3])*(1-x[3,3])+(1-x[2,3])*(1-x[3,3])
# expr_x_line_4 = sum([prod(1 .- prod(x[i,[1,2]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]   # (1-x[1,1]*x[1,2])*(1-x[2,1]*x[2,2])+(1-x[1,1]*x[1,2])*(1-x[3,1]*x[3,2])+(1-x[2,1]*x[2,2])*(1-x[3,1]*x[3,2])
# expr_x_line_5 = sum([prod(1 .- prod(x[i,[1,3]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]   # (1-x[1,1]*x[1,3])*(1-x[2,1]*x[2,3])+(1-x[1,1]*x[1,3])*(1-x[3,1]*x[3,3])+(1-x[2,1]*x[2,3])*(1-x[3,1]*x[3,3])
# expr_x_line_6 = sum([prod(1 .- prod(x[i,[2,3]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1]   # (1-x[1,2]*x[1,3])*(1-x[2,2]*x[2,3])+(1-x[1,2]*x[1,3])*(1-x[3,2]*x[3,3])+(1-x[2,2]*x[2,3])*(1-x[3,2]*x[3,3])
# expr_x_line_7 = sum([prod(1 .- prod(x[i,[1,2,3]],dims=2),dims=1) for i in combinations(1:3,2)])[1,1] # (1-x[1,1]*x[1,2]*x[1,3])*(1-x[2,1]*x[2,2]*x[2,3])+(1-x[1,1]*x[1,2]*x[1,3])*(1-x[3,1]*x[3,2]*x[3,3])+(1-x[2,1]*x[2,2]*x[2,3])*(1-x[3,1]*x[3,2]*x[3,3])
# expr_x = expr_x_line_1 + expr_x_line_2 + expr_x_line_3 - expr_x_line_4 - expr_x_line_5 - expr_x_line_6 + expr_x_line_7
# expand(expr_x)

m,n = 3,40
ki,kj = 1-1,1-1
@eval begin; @syms x[1:$n, 1:$m] δ[1:$n] Δ[1:$n]; end; 
# pfmin = pfmin .- rand(size(pfmin)...)./100

subs_dict = Dict()
for i = 1:n
    # subs_dict[δ[i]] = Float64.(previn[i])
    subs_dict[Δ[i]] = Float64.(pfminneg[i])
    # subs_dict[Δ[i]] = Float64.(1)
    for j = 1:m
        if pfmin[j+kj,i+ki].==1
            # subs_dict[x[i,j]] = Float64.(pfmin[j+kj, i+ki])
        end
        subs_dict[x[i,j]] = Float64.(pfmin[j+kj, i+ki])
        # subs_dict[x[i,j]] = Float64.(1)
    end
end
# subs_dict
pfmin
1 .-pfmin[1+kj:kj+m,1+ki:ki+n]
terms = sum([(-1)^length(j)*prod(1 .- δ[1:n].*(1 .-Δ[1:n].*prod(x[1:n,j],dims=2))[:,1]) for j in powerset(1:m)])
terms_filled = subs(sum([(-1)^length(j)*prod(1 .- δ[1:n].*(1 .-Δ[1:n].*prod(x[1:n,j],dims=2))[:,1]) for j in powerset(1:m)]),subs_dict)
for term in [terms_filled.args...]; println(term,"\n"); end
# expand(subs(terms,subs_dict))
1 .- prod(1 .- rand(size(previn)...)./10 .*previn)
subs_dict
t1 = -(1 - 0.2⋅δ[1]) ⋅(1 - 0.027⋅δ[8])⋅(1 - 0.21⋅δ[19])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
t2 = -(1 - 0.1⋅δ[7]) ⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.1⋅δ[29])
t3 = -(1 - 0.1⋅δ[5]) ⋅(1 - 0.1⋅δ[13])⋅(1 - 0.09⋅δ[17])⋅(1 - 0.09⋅δ[19])⋅(1 - 0.09⋅δ[20])⋅(1 - 0.09⋅δ[22])⋅(1 - 0.09⋅δ[24])⋅(1 - 0.15⋅δ[26])⋅(1 - 0.09⋅δ[28])
t4 = +(1 - 0.2⋅δ[1]) ⋅(1 - 0.1⋅δ[7])⋅(1 - 0.027⋅δ[8])⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.21⋅δ[19])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.1⋅δ[29])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
t5 = +(1 - 0.2⋅δ[1]) ⋅(1 - 0.1⋅δ[5])⋅(1 -  0.027⋅δ[8])⋅(1 - 0.1⋅δ[13])⋅(1 - 0.09⋅δ[17])⋅(1 - 0.2811⋅δ[19])⋅(1 - 0.09⋅δ[20])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.09⋅δ[22])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.09⋅δ[24])⋅(1 - 0.15⋅δ[26])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.09⋅δ[28])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
t6 = +(1 - 0.1⋅δ[5]) ⋅(1 - 0.1⋅δ[7])⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[13])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.09⋅δ[17])⋅(1 - 0.09⋅δ[19])⋅(1 - 0.09⋅δ[20])⋅(1 - 0.09⋅δ[22])⋅(1 - 0.09⋅δ[24])⋅(1 - 0.15⋅δ[26])⋅(1 - 0.09⋅δ[28])⋅(1 - 0.1⋅δ[29])
t7 = -(1 - 0.2⋅δ[1]) ⋅(1 - 0.1⋅δ[5])⋅(1 - 0.1⋅δ[7])⋅(1 - 0.027⋅δ[8]) ⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[13])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.09⋅δ[17])⋅(1 - 0.2811⋅δ[19])⋅(1 - 0.09⋅δ[20])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.09⋅δ[22])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.09⋅δ[24])⋅(1 - 0.15⋅δ[26])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.09⋅δ[28])⋅(1 - 0.1⋅δ[29])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
1 + subs(t1+t2+t3+t4+t5+t6+t7,subs_dict)
subs(t1+t2+t3+t4,subs_dict)
1 .+ subs(t5+t6+t7,subs_dict)
t1 = -(1 - 0.2⋅δ[1]) ⋅(1 - 0.027⋅δ[8])⋅(1 - 0.21⋅δ[19])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
t2 = -(1 - 0.1⋅δ[7]) ⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.1⋅δ[29])
t4 = +(1 - 0.2⋅δ[1]) ⋅(1 - 0.1⋅δ[7])⋅(1 - 0.027⋅δ[8])⋅(1 - 0.1⋅δ[10])⋅(1 - 0.1⋅δ[14])⋅(1 - 0.21⋅δ[19])⋅(1 - 0.081⋅δ[21])⋅(1 - 0.135⋅δ[23])⋅(1 - 0.45⋅δ[27])⋅(1 - 0.1⋅δ[29])⋅(1 - 0.72⋅δ[30])⋅(1 - 0.72⋅δ[35])⋅(1 - 0.9⋅δ[39])
1 .+ subs(t1+t2+t3+t4+t5+t6+t7,subs_dict)
v1 = 1 .- (1 .- pfmin[1,1:40][pfmin[1,1:40].!=1]').*previn[1:40][pfmin[1,1:40].!=1]'
v2 = 1 .- (1 .- pfmin[2,1:40][pfmin[2,1:40].!=1]').*previn[1:40][pfmin[2,1:40].!=1]'
v3 = 1 .- (1 .- pfmin[3,1:40][pfmin[3,1:40].!=1]').*previn[1:40][pfmin[3,1:40].!=1]'
(1 - prod(v1))*(1 - prod(v2))*(1 - prod(v3))
(1 - prod(v1))
(1 - prod(v2))
(1 - prod(v3))

pfmin[1:3,1:40].!=1
pfminneg[1:40]'.!=1
bit = Dict(idx => reduce(|,pfmin[idx,1:40].!=1,dims=1).||(pfminneg[1:40]'.!=1) for idx in powerset(1:3))
terms = Dict(idx => 1 .- prod(pfmin[idx,1:40],dims=1)[bit[idx]]'.*pfminneg[1:40]'[bit[idx]]' for idx in powerset(1:3))
terms_n = Dict(idx => 1 .- prod(pfmin[idx,1:40],dims=1).*pfminneg[1:40]'.*bit[idx] for idx in powerset(1:3))

sommetje = 0
for (i,idx) in enumerate(powerset(1:3))
    # term = (-1)^length(idx).*prod((1 .- terms[idx].*previn[1:40]'[bit[idx]]'))
    term = (-1)^length(idx).*prod((1 .- terms_n[idx].*previn[1:40]'.*bit[idx]))
    sign = (-1)^length(idx)!=1 ? '+' : '-'
    sommetje += term
    println(sign," ",round.((1 .- terms_n[idx].*previn[1:40]'.*bit[idx])[1:20],sigdigits=5))
end
sommetje

terms_filled = subs(sum([(-1)^length(j)*prod(1 .- δ[1:n].*(1 .-Δ[1:n].*prod(x[1:n,j],dims=2))[:,1]) for j in powerset(1:m)]),subs_dict); 
subs_2 = Dict(δ[i]=>previn[i] for i in [2,4,12,15,25,33,34,38]); 
subs_2 = Dict(δ[i]=>previn[i] for i in [2,4,12,15,25,33,34,38,30]); 
subs_2 = Dict(δ[i]=>0 for i in [2,4,12,15,25,33,34,38,1,30,10,13,14,17,20,21,26,22,23,24,35,39]); 
for term in [terms_filled.args...]; println(subs(term,subs_2),"\n"); end
1 .- rand(n)'.*previn[1:n]'
0.999732149544573 * (1-0.2*previn[1]) * 0.1
1-0.1*previn[5]
# 2, 4, 12, 15, 25, 33, 34, 38
x = 1 .- rand(100)./1000
1 .- rand(size(previn)...).*previn


0.999732149544573*(1 - 0.6975*previn[10])*(1 - 0.181*previn[13])*(1 - 0.505*previn[14])*(1 - 0.1719*previn[17])*(1 - 0.09*previn[20])*(1 - 0.099*previn[21])*(1 - 0.54136*previn[26])*(1 - 0.368*previn[29])*(1 - 0.8712*previn[30])*(1 - 0.8712*previn[35])
(1 - previn)
for idx in powerset(1:3)
    println(bit[idx])
end
1 - 0.999732079813255-0.999732076789066





