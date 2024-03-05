include("0. useful/packages.jl")

x = BigInt.([0:100...])
pi_1 =  426880 .*sqrt(BigFloat(10005)) ./ cumsum(factorial.(6 .*x).*(13591409 .+545140134 .*x)./(factorial.(3 .*x).*factorial.(x).^3 .* (-640320).^(3 .*x)))
error_1 = abs.(pi_1 .- pi)./pi
pi_2 = 6 .*cumsum((-1).^x .* (1. ./sqrt(BigFloat(3))).^(2 .*x.+1)./(2 .*x.+1))
error_2 = abs.(pi_2 .- pi)./pi

using SymPy
@syms k,n
expressions = [n]
for r=2:8
    # calculate the coefficients of k^r - (k-1)^r 
    expr = expand(k^r-(k-1)^r)
    coeff = expr.coeff.(k,[0:r-1...])
    # use the coefficients and known expressions to solve the equation 
    solution = (n^r .- dot(coeff[1:end-1],expressions))./coeff[end]
    # add the factorized solution to the list of expressions
    push!(expressions,factor(solution[end]))
end
x = DataFrame(("$i"=>expressions[i] for i=1:8)...)
x[1,"5"]


string.([1:8...])