include("../packages.jl")
include("../quickscore_algorithm.jl")
include("../run_one_time.jl")

previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=18,n_myset=5,language="Julia");
myset = [1,2,4,7,10]; # just like in MATLAB

y_julia = Vector{Any}(undef,5)
y_julia[1] = prod(pfmin[myset, :],dims=1); 
y_julia[2] = 1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev));
y_julia[3] = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)));
y_julia[4] = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
y_julia[5] = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2); 

# Import the same calculated line from Matlab
y_matlab = Vector{Any}(undef,5)
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
matfile = matopen(joinpath(matlab_dir,"variables/precision_comparison/vars.mat"),"r") 
    y_matlab[1] = read(matfile, "a") 
    y_matlab[2] = read(matfile, "b") 
    y_matlab[3] = read(matfile, "c") 
    y_matlab[4] = read(matfile, "d") 
    y_matlab[5] = read(matfile, "e") 
close(matfile) 

headers = ["prod(matrix,dims=1)","()","log()","exp-sum-log()","prod()"]
println("Comparing Julia and Matlab for different single line operations (%.20f)\n")
for i=1:5
    idx = findall(y_matlab[i]-y_julia[i] .!= 0)
    pretty_table(hcat(y_julia[i][idx][1:3],y_matlab[i][idx][1:3])',formatters = ft_printf("%.20f"),header=repeat([headers[i]],3),row_names=["Julia","MATLAB"])
end

# Het probleem is dat prod(matrix,dims=1)[i] anders werkt dan prod(vector[i])
# prod(vector[i]) en matlab[i] zijn hetzelfde
# diagn=8,18,23,30,52,54,63,99 zijn anders
prod(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset, 8])
println("\n###################################### DIFFERENCE WITHOUT BIGFLOAT #############################################"); 
println("# \t Julia prod-matrix, Julia prod-vector \t|||Julia prod-matrix, Matlab prod-matrix"); 
s=0; for i=1:100
    diff_prod_prodvect = prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i])
    diff_prod_matlab = prod(pfmin[myset,:],dims=1)[i] - y_matlab[1][i]
    println("$i \t $diff_prod_prodvect \t\t\t\t\t\t|||$diff_prod_matlab")
    s += diff_prod_matlab
end
println("Total difference: $s")

# Een Float64 heeft 8 bytes(~16 digits), BigFloat heeft 40 bytes(==> ~80 digits) --> Dit zou mogelijk kunnen helpen 
x = BigFloat(10.0^77) 
y = BigFloat(10.0^77)-1
x-y 
# Een Float64 is upgraded naar BigFloat als een BigFloat betrokken is bij de berekening
x = BigFloat(10); println("\nx = BigFloat($x)\tsize = $(sizeof(x)) bytes\t\ttype = $(typeof(x))")
y = 3.0; println("y = $x\t\tsize = $(sizeof(y)) bytes\t\ttype = $(typeof(y))")
z = 0.0; println("z = $x\t\tsize = $(sizeof(z)) bytes\t\ttype = $(typeof(z))")
z = x-y; println("z = x-y = $(z)\t\tsize = $(sizeof(z)) bytes\t\ttype = $(typeof(z)) --> z is upgraded to type = BigFloat")

# De basis van het probleem (verschil in prod-vector en prod-matrix langs as) kan mogelijk worden verholpen door BigFloat te gebruiken
v = [0.7765593409639301 0.6001184280484421 1.0 0.8386605783125091 1.0; 0.8989878262484864 0.8733811867742846 1.0 0.7327799359295427 1.0; 0.647743588618465 0.6036928773336186 0.8387682376832125 0.7178718726114095 0.7783382118777125; 0.6273366589206473 0.5405814766194981 1.0 0.5448994859318418 0.5998139531731465; 0.6217531809701229 0.6539490033957885 0.8675328912964698 0.6171916016781012 0.8755890666242915; 1.0 0.5247440407375039 0.5025303778464986 0.7756979132319572 1.0; 0.6470226323351563 0.5952455922711103 0.7525101796159012 0.789230365665153 1.0; 0.8328467685077987 0.6571857651425702 1.0 0.6676778258545397 0.8328758189728027]
println("\nError abs-value diff w/wo BigFloat for test-matrix v, comparing prod-vector & prod-matrix \nprod(v[:,i]) ==> prod(BigFloat.(v[:,i])) \t\t prod(v,dims=1)[i] ==> prod(BigFloat.(v),dims=1)[i]"); 
for i=axes(v,2)
    if prod(v,dims=1)[i] - prod(v[:,i]) != 0
        prod_diff_BigFloat = prod(BigFloat.(v),dims=1)[i] - prod(BigFloat.(v[:,i]))
        prod_diff_normal = prod(v,dims=1)[i] - prod(v[:,i])
        println("prod(v[:,$i]) \tBigFloat: ",prod_diff_BigFloat,"\tNormal: ",@sprintf("%.0e",abs(prod_diff_normal)))#,"\tRow: ",v[:,i])
    end
end
# Dan ziet het probleem er in de parktijk zo uit
println("\nError abs-value diff w/wo BigFloat for actual problem \nprod(pfmin[myset, i]) ==> prod(BigFloat.(pfmin[myset, i])) \t\t prod(pfmin[myset,:],dims=1)[i] ==> prod(BigFloat.(pfmin[myset,:]),dims=1)[i]"); for i=1:100
    if prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i]) != 0
        prod_diff_BigFloat = prod(BigFloat.(pfmin[myset,:]),dims=1)[i] - prod(BigFloat.(pfmin[myset, i]))
        prod_diff_normal = prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i])
        println("diagn: $i\tBigFloat: ",prod_diff_BigFloat,"\tNormal: ",@sprintf("%.0e",abs(prod_diff_normal)))#,"\tRow: ",pfmin[myset,i])
    end
end

# Lukraak upgraden van prod(matrix,dims=1)[i] via BigFloat helpt echter niet altijd 
BigFloat.(prod(pfmin[myset,:],dims=1)[8]) - BigFloat.(y_matlab[1][8])
prod(BigFloat.(pfmin[myset,:]),dims=1)[8] - BigFloat.(y_matlab[1][8])
# Het upgraden van BigFloat kan ook een nieuw verschil veroorzaken zoals tussen prod(vector[i]) en matlab[i]
BigFloat.(prod(pfmin[myset, 8])) - BigFloat.(y_matlab[1][8])
prod(BigFloat.(pfmin[myset, 8])) - BigFloat.(y_matlab[1][8])
# Echter, door dit verschil zijn prod(matrix,dims=1)[i] en prod(vector[i]) nu wel ineens gelijk 
prod(BigFloat.(pfmin[myset,:]),dims=1)[8] .- prod(BigFloat.(pfmin[myset,8])) 
# Door ze weer terug te brengen via Float64() zijn ze nu hetzelfde net zoals oorspronkelijk
Float64(prod(BigFloat.(pfmin[myset,:]),dims=1)[8]) .- Float64(prod(BigFloat.(pfmin[myset,8])))
# In plaats van het aanpassen voor elk los element kan het via Float64.() ook in een keer gedaan worden 
# Echter, door BigFloat zijn nu matlab[i] en prod(BigFloat(vector[i])) niet meer hetzelfde 
# diagn=29,49,52,59,61,63,72,80,89 zijn anders
Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[8] .- Float64(prod(BigFloat.(pfmin[myset,8])))
println("\n######################################## DIFFERENCE WITH BIGFLOAT ##############################################"); 
println("# \t Julia prod-matrix, Julia prod-vector \t\t|||Julia prod-matrix, Matlab prod-matrix"); 
s=0; for i=1:100
    diff_prod_prodvect = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - Float64(prod(BigFloat.(pfmin[myset,i])))
    diff_prod_matlab = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - y_matlab[1][i]
    println("$i \t $diff_prod_prodvect \t\t\t\t\t\t|||$diff_prod_matlab")
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
prod1(pfmin[myset,:],dims=1)[8] - y_matlab[1][8]
prod1(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset,:],dims=1)[8]
prod(pfmin[myset,:],dims=1)[8] - y_matlab[1][8]
println("\n######################################## DIFFERENCE WITH 'prod1' ##############################################"); 
println("# \t Julia prod1-matrix, Matlab prod-matrix \t|||Julia prod1-matrix, Julia prod-matrix"); 
for i=1:100
    diff_prod1_prod = prod(pfmin[myset,:],dims=1)[i] - prod1(pfmin[myset,:],dims=1)[i]
    diff_prod1_matlab = prod1(pfmin[myset,:],dims=1)[i] - y_matlab[1][i]
    println("$i \t $diff_prod1_matlab \t\t\t\t\t\t $diff_prod1_prod")
end

# Vergelijking van tijden prod1 en prod geeft wel dat prod1 3 keer langzamer is dan prod 
# @benchmark prod1(pfmin[myset,:],dims=1)
# @benchmark prod(pfmin[myset,:],dims=1)

# Nu kan ik kijken naar de verschillende implementatie regels a-e (ik heb f toegevoegd nog)
y_julia_prod1 = Vector{Any}(undef,6);
y_julia_prod1[1] = prod1(pfmin[myset, :],dims=1);
y_julia_prod1[2] = 1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev));
y_julia_prod1[3] = log.(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)));
y_julia_prod1[4] = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
y_julia_prod1[5] = ((-1)^length(myset)) .* prod1(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2);
y_julia_prod1[6] = ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2);
# Het is duidelijk te zien dat de stap waarbij de log genomen wordt verschillend is nog nu 
println("\na-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[1] == y_matlab[1]) \t-->\t a = prod1(pfmin[myset, :],dims=1)")
println("b-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[2] == y_matlab[2]) \t-->\t b = 1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))")
println("c-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[3] == y_matlab[3]) \t-->\t c = log.(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)))")
println("d-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[4] == y_matlab[4]) \t-->\t d = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2))")
println("e-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[5] == y_matlab[5]) \t-->\t e = ((-1)^length(myset)) .* prod1(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2)")
println("f-matrix for prod1 is same for Julia and Matlab:\t$(y_julia_prod1[6] == y_matlab[5]) \t-->\t f = ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2)")
pretty_table(hcat(y_julia[2][1,62],y_matlab[2][1,62])',formatters = ft_printf("%.50f"),header=["b = 1e-50 .+ (a .* prevminneg .+ (1 .- prev))  at [1,62]"],row_names=["Julia","MATLAB"],alignment=:l)
pretty_table(hcat(y_julia[3][1,62],y_matlab[3][1,62])',formatters = ft_printf("%.50f"),header=["c = log.(b)                                    at [1,62]"],row_names=["Julia","MATLAB"],alignment=:l)

# Het blijkt dat log(getal) niet altijd hetzelfde is voor een simpel getal in Julia en MATLAB (gekopieerd uit MATLAB-folder)
# Dit komt zeer waarschijnlijk door de afronding na 16 digits
getal = 0.97150353070211348
x_test = log(getal)
x_test_matlab = -0.02891037589217936021
pretty_table(hcat(vcat(getal, getal),vcat(x_test,x_test_matlab)),formatters = ft_printf("%.20f"),header=["getal","log(getal)"],row_names=["Julia","Matlab"])
println("Conclusion: log(getal) starts to differ in Julia and Matlab after 16 digits, because it is Float64")

# We laten log.() voor nu even en kijken naar prod() + prod1 zoals in f --> verschil QS Matlab en Julia QS zou dus 0 moeten zijn
matfile = matopen(joinpath(matlab_dir,"variables/random_m_1_22/m_18.mat"),"r")
    posterior_matlab_exp_sum_log = read(matfile, "posterior")'
    pfplus_matlab_exp_sum_log = read(matfile,"pfplus")
close(matfile)
previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=18,n_myset=3,language="Julia");
pfplus_julia_prod, P_joint_julia_prod, posterior_julia_prod,dt_julia_prod = quickscore(previn, pfmin, pfminneg,"prod prod1");  
pfplus_julia_exp_sum_log, P_joint_julia_exp_sum_log, posterior_julia_exp_sum_log,dt_julia_exp_sum_log = quickscore(previn, pfmin, pfminneg,"exp-sum-log prod1");  

# Een vergelijking binnen (en tussen) Matlab en Julia voor prod en log laat zien dat prod tussen Julia en Matlab het kleinste verschil is
println("\n################################### Posterior comparison using prod1(pfmin) ###################################")
pretty_table(hcat((posterior_julia_prod - posterior_julia_exp_sum_log),(posterior_matlab - posterior_julia_prod),(posterior_matlab - posterior_julia_exp_sum_log)),header=["Julia prod  -  Julia exp-sum-log","MATLAB prod  -  Julia prod","MATLAB exp-sum-log  -  Julia exp-sum-log"])
println("\nConclusion 1: After 16 digits precision stops and in the end this causes precision problems after lots of calculations like prod/exp-sum-log and summing. BigFloat can help here. ")
println("\nConclusion 2: prod-matrix and prod1-matrix have just different values after the cutoff of 16 digits. It is a 'coincidence' that prod1-matrix is the same as MATLAB's prod-matrix.")






