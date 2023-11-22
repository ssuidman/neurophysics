include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")

function quickscore_Fortran(previn,pfmin,pfminneg;method="prod")
    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $m (postive tests)")
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    pfplus = zeros(Float64,n+1,1)  # will be used for P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1) (note F- will be absorbed below) 

    dt = @elapsed begin
        for i in ProgressBar(0:(2^m-1)) 
            v = digits(i,base=2,pad=m) 
            myset = findall(v.==1) 
            # pfplus = pfplus + ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
            # pfplus = pfplus + ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
            pfplus = pfplus + ((-1)^length(myset)) .* prod(1e-50 .+ (prod(BigFloat.(pfmin[myset, :]),dims=1) .* prevminneg .+ (1 .- prev)),dims=2); 
            # pfplus = pfplus + ((-1)^length(myset)) .* prod(1e-50 .+ (BigFloat.(pfmin[[1], :]) .* prevminneg .+ (1 .- prev)),dims=2); 
        end
    end
    println("\nRunning time: $dt")

    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    return P_joint, posterior, dt 
end

P_joint, posterior, dt = quickscore_Fortran(previn,pfmin,pfminneg)



previn = Matrix(DataFrame(CSV.File("variables/cpp/patient404_case_1_previn.csv",header=false)))[:,1]
pfmin = Matrix(DataFrame(CSV.File("variables/cpp/patient404_case_1_pfmin.csv",header=false)))
pfminneg = Matrix(DataFrame(CSV.File("variables/cpp/patient404_case_1_pfminneg.csv",header=false)))[:,1]
P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"prod BF")




precision(Float128)
x = rand(100,101)
y = BigFloat.(x,precision=53)
z = Float128.(rand(100,101))

@btime prod(x,dims=1)
@btime prod(y,dims=1)
@btime prod(z,dims=1)

include("../packages.jl")



