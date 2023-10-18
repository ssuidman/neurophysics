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
# pfplus,posterior,dt = quickscore(previn, pfmin, pfminneg,false); 

# Set some myset values
myset = [1,2,4,7,10]

# Reken een aantal termen uit in Julia gesorteerd van complex naar simpel
a_j = prod(pfmin[myset, :],dims=1)
b_j = 1e-50 .+ (a_j .* prevminneg .+ (1 .- prev))
c_j = log.(b_j)
d_j = ((-1)^length(myset)) .* exp.(sum(c_j, dims=2))
e_j = ((-1)^length(myset)) .* prod(b_j, dims=2)

# Import the same calculated line from Matlab
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables_julia/test.mat"),"r") 
a_m = read(matfile, "a") 
b_m = read(matfile, "b") 
c_m = read(matfile, "c") 
d_m = read(matfile, "d") 
e_m = read(matfile, "e") 
close(matfile) 

# There is a difference between Julia and Matlab, even for very basic function such as in a_j and a_m --> prod(matrix,dims=1)
a_i = findall(a_m-a_j .!= 0)
pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), a_j[a_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), a_m[a_i][1:3])),header=String["Julia prod(matrix,dims=1)","Matlab prod(matrix,dims=1)"])
b_i = findall(b_m-b_j .!= 0)
pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), b_j[b_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), b_m[b_i][1:3])),header=String["Julia ()","Matlab ()"])
c_i = findall(c_m-c_j .!= 0)
pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), c_j[c_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), c_m[c_i][1:3])),header=String["Julia log()","Matlab log()"])
d_i = findall(d_m-d_j .!= 0)
pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), d_j[d_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), d_m[d_i][1:3])),header=String["Julia exp-sum-log()","Matlab exp-sum-log()"])
e_i = findall(e_m-e_j .!= 0)
pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), e_j[e_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), e_m[e_i][1:3])),header=String["Julia prod()","Matlab prod()"])

# Het probleem is dat prod(matrix,dims=1)[i] anders werkt dan prod(vector[i])
# prod(vector[i]) en matlab[i] zijn hetzelfde
# diagn=8,18,23,30,52,54,63,99 zijn anders
prod(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset, 8])
prod(pfmin[myset,:],dims=1)[8] - a_m[8]
prod(pfmin[myset, 8]) - a_m[8]
s=0; for i=1:100
    diff_prod_prodvect = prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i])
    diff_prod_matlab = prod(pfmin[myset,:],dims=1)[i] - a_m[i]
    println("$i \t $diff_prod_prodvect \t $diff_prod_matlab")
    s += diff_prod_matlab
end
println("Total difference: $s")

# Een Float64 heeft 8 bytes(~16 digits), BigFloat heeft 40 bytes(==> ~80 digits) --> Dit zou mogelijk kunnen helpen 
x = BigFloat(10.0^77) 
y = BigFloat(10.0^77)-1
x-y 

# Een Float64 is upgraded naar BigFloat als een BigFloat betrokken is bij de berekening
x = BigFloat(10); print("x = BigFloat($x)\tsize = $(sizeof(x)) bytes\t\ttype = $(typeof(x))")
y = 3.0; print("y = $x\t\tsize = $(sizeof(y)) bytes\t\ttype = $(typeof(y))")
z = 0.0; print("z = $x\t\tsize = $(sizeof(z)) bytes\t\ttype = $(typeof(z))")
z = x-y; print("z = x-y = $(z)\t\tsize = $(sizeof(z)) bytes\t\ttype = $(typeof(z)) --> z is upgraded to type = BigFloat\n")

# De basis van het probleem kan mogelijk worden verholpen door BigFloat te gebruiken
v = [0.7765593409639301 0.6001184280484421 1.0 0.8386605783125091 1.0; 0.8989878262484864 0.8733811867742846 1.0 0.7327799359295427 1.0; 0.647743588618465 0.6036928773336186 0.8387682376832125 0.7178718726114095 0.7783382118777125; 0.6273366589206473 0.5405814766194981 1.0 0.5448994859318418 0.5998139531731465; 0.6217531809701229 0.6539490033957885 0.8675328912964698 0.6171916016781012 0.8755890666242915; 1.0 0.5247440407375039 0.5025303778464986 0.7756979132319572 1.0; 0.6470226323351563 0.5952455922711103 0.7525101796159012 0.789230365665153 1.0; 0.8328467685077987 0.6571857651425702 1.0 0.6676778258545397 0.8328758189728027]
println("Error absolute value difference with/without BigFloat:"); for i=axes(v,2)
    if prod(v,dims=1)[i] - prod(v[:,i]) != 0
        prod_diff_BigFloat = prod(BigFloat.(v),dims=1)[i] - prod(BigFloat.(v[:,i]))
        prod_diff_normal = prod(v,dims=1)[i] - prod(v[:,i])
        println("prod(v[:,$i]) \tBigFloat: ",prod_diff_BigFloat,"\tNormal: ",@sprintf("%.0e",abs(prod_diff_normal)),"\tRow: ",v[:,i])
    end
end
# Dan ziet het probleem er in de parktijk zo uit
println("Error absolute value difference with/without BigFloat:"); for i=1:100
    if prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i]) != 0
        prod_diff_BigFloat = prod(BigFloat.(pfmin[myset,:]),dims=1)[i] - prod(BigFloat.(pfmin[myset, i]))
        prod_diff_normal = prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i])
        println("diagn: $i\tBigFloat: ",prod_diff_BigFloat,"\tNormal: ",@sprintf("%.0e",abs(prod_diff_normal)),"\tRow: ",pfmin[myset,i])
    end
end

# Lukraak upgraden van prod(matrix,dims=1)[i] via BigFloat helpt echter niet altijd 
BigFloat.(prod(pfmin[myset,:],dims=1)[8]) - BigFloat.(a_m[8])
prod(BigFloat.(pfmin[myset,:]),dims=1)[8] - BigFloat.(a_m[8])
# Het upgraden van BigFloat kan ook een nieuw verschil veroorzaken zoals tussen prod(vector[i]) en matlab[i]
BigFloat.(prod(pfmin[myset, 8])) - BigFloat.(a_m[8])
prod(BigFloat.(pfmin[myset, 8])) - BigFloat.(a_m[8])
# Echter, door dit verschil zijn prod(matrix,dims=1)[i] en prod(vector[i]) nu wel ineens gelijk 
prod(BigFloat.(pfmin[myset,:]),dims=1)[8] .- prod(BigFloat.(pfmin[myset,8])) 
for i=1:100; println("$i \t",prod(BigFloat.(pfmin[myset,:]),dims=1)[i] .- prod(BigFloat.(pfmin[myset,i]))); end
# Door ze weer terug te brengen via Float64() zijn ze nu gelijk 
Float64(prod(BigFloat.(pfmin[myset,:]),dims=1)[8]) .- Float64(prod(BigFloat.(pfmin[myset,8])))
for i=1:100
    diff_prod_prodvect = Float64(prod(BigFloat.(pfmin[myset,:]),dims=1)[i]) - Float64(prod(BigFloat.(pfmin[myset,i])))
    println("$i \t",diff_prod_prodvect)
end
# In plaats van het aanpassen voor elk los element kan het via Float64.() ook in een keer gedaan worden 
# Echter, door BigFloat zijn nu matlab[i] en prod(BigFloat(vector[i])) niet meer hetzelfde 
# diagn=29,49,52,59,61,63,72,80,89 zijn anders
Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[8] .- Float64(prod(BigFloat.(pfmin[myset,8])))
s=0; for i=1:100
    diff_prod_prodvect = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - Float64(prod(BigFloat.(pfmin[myset,i])))
    diff_prod_matlab = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - a_m[i]
    println("$i \t $diff_prod_prodvect \t $diff_prod_matlab")
    s += diff_prod_matlab
end
println("Total difference: $s")

# We kunnen een functie maken prod1 die voor een matrix hetzelfde werkt als prod 
function prod1(x;dims=1)
    if dims==1
        result = zeros(1,size(x,2))
        for i=axes(x,2)
            result[1,i] = prod(x[:,i])
        end
    elseif dims==2
        result = zeros(size(x,1),1)
        for i=axes(x,1)
            result[i,1] = prod(x[i,:])
        end
    end    
    return result
end

# Een precisievergelijking laat zien dat prod1 hetzelfde geeft als matlab, maar niet als prod 
prod1(pfmin[myset,:],dims=1)[8] - a_m[8]
prod1(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset,:],dims=1)[8]
prod(pfmin[myset,:],dims=1)[8] - a_m[8]
for i=1:100
    diff_prod1_matlab = prod1(pfmin[myset,:],dims=1)[i] - a_m[i]
    diff_prod1_prod = prod(pfmin[myset,:],dims=1)[i] - prod1(pfmin[myset,:],dims=1)[i]
    println("$i \t $diff_prod1_matlab \t $diff_prod1_prod")
end

# Vergelijking van tijden prod1 en prod geeft wel dat prod1 3 keer langzamer is dan prod 
@benchmark prod1(pfmin[myset,:],dims=1)
@benchmark prod(pfmin[myset,:],dims=1)

# Nu kan ik kijken naar de verschillende implementatie regels a-e (ik heb f toegevoegd nog)
a1_j = prod1(pfmin[myset, :],dims=1)
b1_j = 1e-50 .+ (a1_j .* prevminneg .+ (1 .- prev))
c1_j = log.(b1_j)
d1_j = ((-1)^length(myset)) .* exp.(sum(c1_j, dims=2))
e1_j = ((-1)^length(myset)) .* prod1(b1_j, dims=2)
f1_j = ((-1)^length(myset)) .* prod(b1_j, dims=2)
# Het is duidelijk te zien dat de stap waarbij de log genomen wordt verschillend is nog nu 
println("a-matrix for prod1 is same for Julia and Matlab:\t$(a1_j == a_m)")
println("b-matrix for prod1 is same for Julia and Matlab:\t$(b1_j == b_m)")
println("c-matrix for prod1 is same for Julia and Matlab:\t$(c1_j == c_m)")
println("d-matrix for prod1 is same for Julia and Matlab:\t$(d1_j == d_m)")
println("e-matrix for prod1 is same for Julia and Matlab:\t$(e1_j == e_m)")
println("f-matrix for prod1 is same for Julia and Matlab:\t$(f1_j == e_m)")
pretty_table(vcat(Printf.format.(Ref(Printf.Format("%.70f")), b1_j[1,62]),Printf.format.(Ref(Printf.Format("%.70f")), b_m[1,62])),header=String["()       Boven: Julia        Onder: Matlab"])
pretty_table(vcat(Printf.format.(Ref(Printf.Format("%.70f")), c1_j[1,62]),Printf.format.(Ref(Printf.Format("%.70f")), c_m[1,62])),header=String["log()       Boven: Julia        Onder: Matlab"])

# Het blijkt dat log(getal) niet altijd hetzelfde is voor een simpel getal in Julia en MATLAB (gekopieerd uit MATLAB-folder)
# Dit komt zeer waarschijnlijk door de afronding na 16 digits
getal = 0.97150353070211348
x_test = log(getal)
x_test_matlab = -0.02891037589217936021
print("Julia: \t log(x_test)= \n",Printf.format.(Ref(Printf.Format("%.20f")), x_test))
print("MATLAB:\t log(x_test)= \n",Printf.format.(Ref(Printf.Format("%.20f")), x_test_matlab))


