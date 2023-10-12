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



# Verschil tussen matlab en julia in:   "prod(pfmin[myset, :],dims=1)"
prod(BigFloat.(pfmin[myset, :]),dims=1)[8] - BigFloat.(a_m[8])
BigFloat.(prod(pfmin[myset, :],dims=1)[8]) - BigFloat.(a_m[8])
BigFloat.(prod(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset, 8]))
prod(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset, 8])
# Deze is wel hetzelfde (verschil=0), zoals het zou moeten 
prod(pfmin[myset, 8]) - a_m[8] 
BigFloat.(prod(pfmin[myset, 8])) - a_m[8]

x = BigFloat(10.0^77) # maximale macht is 77 waarbij het nog steeds werkt 
y = BigFloat(10.0^77)-1 # zelfde maar dan min 1
x-y # kijken naar het verschil

x = BigFloat(10)
sizeof(x)
typeof(x)
y = 3.0
sizeof(y)
typeof(y)
z = 0.0
sizeof(z)
typeof(z)
z = x-y
sizeof(z)
typeof(z)




prod(BigFloat.(pfmin[myset,:]),dims=1)[8] - BigFloat.(a_m[8])

for k=1:1000
    Random.seed!(1234) # first run this 
    x = rand(5,100) # then run this 
    s = 0
    for i=1:100
        s += prod(x,dims=1)[8]-prod(x[:,8])
    end
end

for (j,i) in enumerate(1:100)
    pfmin[myset,:][:,i]
    if prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i]) != 0
        println(j,pfmin[myset,ii],prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i]))
    end
end

Random.seed!(1234) # first run this 
v_rand = rand(5,8) 

v = hcat([0.7765593409639301, 0.6001184280484421, 1.0, 0.8386605783125091, 1.0], [0.8989878262484864, 0.8733811867742846, 1.0, 0.7327799359295427, 1.0], [0.647743588618465, 0.6036928773336186, 0.8387682376832125, 0.7178718726114095, 0.7783382118777125], [0.6273366589206473, 0.5405814766194981, 1.0, 0.5448994859318418, 0.5998139531731465], [0.6217531809701229, 0.6539490033957885, 0.8675328912964698, 0.6171916016781012, 0.8755890666242915], [1.0, 0.5247440407375039, 0.5025303778464986, 0.7756979132319572, 1.0], [0.6470226323351563, 0.5952455922711103, 0.7525101796159012, 0.789230365665153, 1.0], [0.8328467685077987, 0.6571857651425702, 1.0, 0.6676778258545397, 0.8328758189728027]) 
prod(v,dims=1) .- prod(v[:,1]) 
prod(BigFloat.(v),dims=1) .- prod(BigFloat.(v)[:,1]) 

