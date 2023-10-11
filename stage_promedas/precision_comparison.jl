include("packages.jl")
include("quickscore_preparation.jl")
include("quickscore_algorithm.jl")

# Load random variables
file = jldopen("variables/m_18.jld", "r")
previn = file["previn"]
pfminneg = file["pfminneg"]
pfmin = file["pfmin"]
close(file)

# Set some variables based on the loaded variables
m,n = size(pfmin)
prev = repeat(previn',inner=(n+1,1)); for i in 1:n prev[i+1,i]=1 end 
prevend = vcat(1,previn) # prevend are the prevalences 1, P(d_1=1),... P(d_n=1) and is needed to multiply at the end
prevminneg = prev.*pfminneg'

# Run the quickscore algorithm
pfplus,posterior,dt = quickscore(previn, pfmin, pfminneg,false); 

# Set some myset values
myset = [1,2,4,7,10]

# Calculate a line in Julia
x_j = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))
y_j = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))
z_j = 1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))
# z_j = 1e-50 .+ (reshape(previn,1,100) .* prevminneg .+ (1 .- prev))
a_j = prod(pfmin[myset, :],dims=1)

# Import the same calculated line from Matlab
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/test.mat"),"r") 
x_m = read(matfile, "x") 
y_m = read(matfile, "y") 
z_m = read(matfile, "z") 
a_m = read(matfile, "a") 
close(matfile) 

# printed_y = Printf.format.(Ref(Printf.Format("%.70f")), y_julia)
# printed_y[1:10,1:10]

x_diff = x_m-x_j
y_diff = y_m-y_j
z_diff = z_m-z_j
a_diff = a_m-a_j

x_i = findall(x_diff .!= 0)
Printf.format.(Ref(Printf.Format("%.70f")), x_m[x_i][1:3])
Printf.format.(Ref(Printf.Format("%.70f")), x_j[x_i][1:3])

y_i = findall(y_diff .!= 0)
Printf.format.(Ref(Printf.Format("%.70f")), y_m[y_i][1:3])
Printf.format.(Ref(Printf.Format("%.70f")), y_j[y_i][1:3])

z_i = findall(z_diff .!= 0)
Printf.format.(Ref(Printf.Format("%.70f")), z_m[z_i][1:3])
Printf.format.(Ref(Printf.Format("%.70f")), z_j[z_i][1:3])

a_i = findall(a_diff .!= 0)
Printf.format.(Ref(Printf.Format("%.70f")), a_m[a_i][1:3])
Printf.format.(Ref(Printf.Format("%.70f")), a_j[a_i][1:3])

prod(pfmin[myset,8])

# Vergelijkingen met BigFloat
prod(BigFloat.(pfmin[myset, :]),dims=1)[8]
BigFloat.(prod(pfmin[myset, :],dims=1)[8])
BigFloat.(prod(pfmin[myset, 8]))
BigFloat.(a_j[8])
BigFloat.(a_m[8])
# Deze is wel hetzelfde (verschil=0), zoals het zou moeten 
prod(pfmin[myset, 8]) - a_m[8]

