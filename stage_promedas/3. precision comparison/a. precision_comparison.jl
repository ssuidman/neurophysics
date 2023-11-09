include("../packages.jl")
include("../quickscore_algorithm.jl")
include("../run_one_time.jl")

previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=18,n_myset=5,language="Julia");
myset = [1,2,4,7,10]; # just like in MATLAB
# pfplus_, P_joint_, posterior_,dt_ = quickscore(previn, pfmin, pfminneg,"prod");  


y_julia = Vector{Any}(undef,5)
y_julia[1] = prod(pfmin[myset, :],dims=1); 
y_julia[2] = 1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev));
y_julia[3] = log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)));
y_julia[4] = ((-1)^length(myset)) .* exp.(sum(log.(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev))), dims=2));
y_julia[5] = ((-1)^length(myset)) .* prod(1e-50 .+ (prod(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)), dims=2); 


a_j = prod(pfmin[myset, :],dims=1)
b_j = 1e-50 .+ (a_j .* prevminneg .+ (1 .- prev))
c_j = log.(b_j)
d_j = ((-1)^length(myset)) .* exp.(sum(c_j, dims=2))
e_j = ((-1)^length(myset)) .* prod(b_j, dims=2)

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

i = 1

println("Comparing Julia and Matlab for different single line operations\n")
headers = [["Julia prod(matrix,dims=1)","Matlab prod(matrix,dims=1)"],["Julia ()","Matlab ()"],["Julia log()","Matlab log()"],["Julia exp-sum-log()","Matlab exp-sum-log()"],["Julia prod()","Matlab prod()"]]
for i=1:5
    idx = findall(y_matlab[i]-y_julia[i] .!= 0)
    pretty_table(hcat(y_julia[i][idx][1:3],y_matlab[i][idx][1:3])',formatters = ft_printf("%.20f"),row_names=headers[i])
end

# # There is a difference between Julia and Matlab, even for very basic function such as in a_j and a_m --> prod(matrix,dims=1)
# a_i = findall(a_m-a_j .!= 0)
# pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), a_j[a_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), a_m[a_i][1:3])),header=String["Julia prod(matrix,dims=1)","Matlab prod(matrix,dims=1)"])
# b_i = findall(b_m-b_j .!= 0)
# pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), b_j[b_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), b_m[b_i][1:3])),header=String["Julia ()","Matlab ()"])
# c_i = findall(c_m-c_j .!= 0)
# pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), c_j[c_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), c_m[c_i][1:3])),header=String["Julia log()","Matlab log()"])
# d_i = findall(d_m-d_j .!= 0)
# pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), d_j[d_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), d_m[d_i][1:3])),header=String["Julia exp-sum-log()","Matlab exp-sum-log()"])
# e_i = findall(e_m-e_j .!= 0)
# pretty_table(hcat(Printf.format.(Ref(Printf.Format("%.70f")), e_j[e_i][1:3]),Printf.format.(Ref(Printf.Format("%.70f")), e_m[e_i][1:3])),header=String["Julia prod()","Matlab prod()"])

# Het probleem is dat prod(matrix,dims=1)[i] anders werkt dan prod(vector[i])
# prod(vector[i]) en matlab[i] zijn hetzelfde
# diagn=8,18,23,30,52,54,63,99 zijn anders
prod(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset, 8])
prod(pfmin[myset,:],dims=1)[8] - a_m[8]
prod(pfmin[myset, 8]) - a_m[8]
println("\n###################################### DIFFERENCE WITHOUT BIGFLOAT #############################################"); 
println("# \t Julia prod(matrix,dims=1), Julia prod(vector) \t|||Julia prod(matrix,dims=1), Matlab prod(matrix,dims=1)"); 
s=0; for i=1:100
    diff_prod_prodvect = prod(pfmin[myset,:],dims=1)[i] - prod(pfmin[myset, i])
    diff_prod_matlab = prod(pfmin[myset,:],dims=1)[i] - a_m[i]
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

# De basis van het probleem kan mogelijk worden verholpen door BigFloat te gebruiken
v = [0.7765593409639301 0.6001184280484421 1.0 0.8386605783125091 1.0; 0.8989878262484864 0.8733811867742846 1.0 0.7327799359295427 1.0; 0.647743588618465 0.6036928773336186 0.8387682376832125 0.7178718726114095 0.7783382118777125; 0.6273366589206473 0.5405814766194981 1.0 0.5448994859318418 0.5998139531731465; 0.6217531809701229 0.6539490033957885 0.8675328912964698 0.6171916016781012 0.8755890666242915; 1.0 0.5247440407375039 0.5025303778464986 0.7756979132319572 1.0; 0.6470226323351563 0.5952455922711103 0.7525101796159012 0.789230365665153 1.0; 0.8328467685077987 0.6571857651425702 1.0 0.6676778258545397 0.8328758189728027]
println("\nError abs-value diff w/wo BigFloat for test-matrix v \nprod(v[:,i]) ==> prod(BigFloat.(v[:,i])) \t\t prod(v,dims=1)[i] ==> prod(BigFloat.(v),dims=1)[i]"); for i=axes(v,2)
    if prod(v,dims=1)[i] - prod(v[:,i]) != 0
        prod_diff_BigFloat = prod(BigFloat.(v),dims=1)[i] - prod(BigFloat.(v[:,i]))
        prod_diff_normal = prod(v,dims=1)[i] - prod(v[:,i])
        println("prod(v[:,$i]) \tBigFloat: ",prod_diff_BigFloat,"\tNormal: ",@sprintf("%.0e",abs(prod_diff_normal)),"\tRow: ",v[:,i])
    end
end
# Dan ziet het probleem er in de parktijk zo uit
println("\nError abs-value diff w/wo BigFloat for actual problem \nprod(pfmin[myset, i]) ==> prod(BigFloat.(pfmin[myset, i])) \t\t prod(pfmin[myset,:],dims=1)[i] ==> prod(BigFloat.(pfmin[myset,:]),dims=1)[i]"); for i=1:100
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
# Door ze weer terug te brengen via Float64() zijn ze nu hetzelfde net zoals oorspronkelijk
Float64(prod(BigFloat.(pfmin[myset,:]),dims=1)[8]) .- Float64(prod(BigFloat.(pfmin[myset,8])))
# In plaats van het aanpassen voor elk los element kan het via Float64.() ook in een keer gedaan worden 
# Echter, door BigFloat zijn nu matlab[i] en prod(BigFloat(vector[i])) niet meer hetzelfde 
# diagn=29,49,52,59,61,63,72,80,89 zijn anders
Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[8] .- Float64(prod(BigFloat.(pfmin[myset,8])))
println("\n######################################## DIFFERENCE WITH BIGFLOAT ##############################################"); 
println("# \t Julia prod(matrix,dims=1), Julia prod(vector) \t|||Julia prod(matrix,dims=1), Matlab prod(matrix,dims=1)"); 
s=0; for i=1:100
    diff_prod_prodvect = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - Float64(prod(BigFloat.(pfmin[myset,i])))
    diff_prod_matlab = Float64.(prod(BigFloat.(pfmin[myset,:]),dims=1))[i] - a_m[i]
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
prod1(pfmin[myset,:],dims=1)[8] - a_m[8]
prod1(pfmin[myset,:],dims=1)[8] - prod(pfmin[myset,:],dims=1)[8]
prod(pfmin[myset,:],dims=1)[8] - a_m[8]
println("\n######################################## DIFFERENCE WITH 'prod1' ##############################################"); 
println("# \t Julia prod1(pfmin), Matlab prod(pfmin) \t|||Julia prod1(pfmin), Julia prod(pfmin)"); 
for i=1:100
    diff_prod1_prod = prod(pfmin[myset,:],dims=1)[i] - prod1(pfmin[myset,:],dims=1)[i]
    diff_prod1_matlab = prod1(pfmin[myset,:],dims=1)[i] - a_m[i]
    println("$i \t $diff_prod1_matlab \t\t\t\t\t\t $diff_prod1_prod")
end

# Vergelijking van tijden prod1 en prod geeft wel dat prod1 3 keer langzamer is dan prod 
# @benchmark prod1(pfmin[myset,:],dims=1)
# @benchmark prod(pfmin[myset,:],dims=1)

# Nu kan ik kijken naar de verschillende implementatie regels a-e (ik heb f toegevoegd nog)
a1_j = prod1(pfmin[myset, :],dims=1)
b1_j = 1e-50 .+ (a1_j .* prevminneg .+ (1 .- prev))
c1_j = log.(b1_j)
d1_j = ((-1)^length(myset)) .* exp.(sum(c1_j, dims=2))
e1_j = ((-1)^length(myset)) .* prod1(b1_j, dims=2)
f1_j = ((-1)^length(myset)) .* prod(b1_j, dims=2)
# Het is duidelijk te zien dat de stap waarbij de log genomen wordt verschillend is nog nu 
println("\na-matrix for prod1 is same for Julia and Matlab:\t$(a1_j == a_m) \t-->\t a = prod1(pfmin[myset, :],dims=1)")
println("b-matrix for prod1 is same for Julia and Matlab:\t$(b1_j == b_m) \t-->\t b = 1e-50 .+ (a .* prevminneg .+ (1 .- prev))")
println("c-matrix for prod1 is same for Julia and Matlab:\t$(c1_j == c_m) \t-->\t c = log.(b)")
println("d-matrix for prod1 is same for Julia and Matlab:\t$(d1_j == d_m) \t-->\t d = ((-1)^length(myset)) .* exp.(sum(c, dims=2))")
println("e-matrix for prod1 is same for Julia and Matlab:\t$(e1_j == e_m) \t-->\t e = ((-1)^length(myset)) .* prod1(b, dims=2)")
println("f-matrix for prod1 is same for Julia and Matlab:\t$(f1_j == e_m) \t-->\t f = ((-1)^length(myset)) .* prod(b, dims=2)")
pretty_table(vcat(Printf.format.(Ref(Printf.Format("%.70f")), b1_j[1,62]),Printf.format.(Ref(Printf.Format("%.70f")), b_m[1,62])),header=String["b = 1e-50 .+ (a .* prevminneg .+ (1 .- prev))         at [1,62]"],row_names=["Julia","Matlab"])
pretty_table(vcat(Printf.format.(Ref(Printf.Format("%.70f")), c1_j[1,62]),Printf.format.(Ref(Printf.Format("%.70f")), c_m[1,62])),header=String["c = log.(b)                                           at [1,62]"],row_names=["Julia","Matlab"])
# Het blijkt dat log(getal) niet altijd hetzelfde is voor een simpel getal in Julia en MATLAB (gekopieerd uit MATLAB-folder)
# Dit komt zeer waarschijnlijk door de afronding na 16 digits
getal = 0.97150353070211348
x_test = log(getal)
x_test_matlab = -0.02891037589217936021
pretty_table(hcat(vcat(Printf.format.(Ref(Printf.Format("%.20f")), getal),Printf.format.(Ref(Printf.Format("%.20f")), getal)),vcat(Printf.format.(Ref(Printf.Format("%.20f")), x_test),Printf.format.(Ref(Printf.Format("%.20f")), x_test_matlab))),header=String["getal","log(getal)"],row_names=["Julia","Matlab"])
println("Conclusion: log(getal) starts to differ in Julia and Matlab after exactly 16 digits, because it is Float64")

# We laten log.() voor nu even en kijken naar prod() + prod1 zoals in f --> verschil QS Matlab en Julia QS zou dus 0 moeten zijn
matfile = matopen(joinpath(matlab_dir,"variables_julia/m_18.mat"),"r")
    posterior_prod_matlab = read(matfile, "posterior_prod")'
    posterior_log_matlab = read(matfile, "posterior_log")'
    pfplus_matlab = read(matfile,"pfplus")
    pfplus_matrix_matlab = read(matfile,"pfplus_matrix")
close(matfile)
jldopen("variables/m_18.jld", "r") do file
    global posterior_prod_julia = file["posterior_prod"]
    global posterior_log_julia = file["posterior_log"] 
end
# Een vergelijking binnen (en tussen) Matlab en Julia voor prod en log laat zien dat prod tussen Julia en Matlab het kleinste verschil is
println("\n################################### Posterior comparison (using prod1[pfmin]) ###################################")
pretty_table(hcat((posterior_prod_julia - posterior_log_julia),(posterior_prod_matlab-posterior_log_matlab),(posterior_prod_matlab - posterior_prod_julia),(posterior_log_matlab - posterior_log_julia)),header=["Julia prod  -  Julia log","MATLAB prod  -  MATLAB log","MATLAB prod  -  Julia prod","MATLAB log  -  Julia log"])

# We proberen nu pfplus_matrix te maken in Julia en Matlab om te kijken waar het nog anders is. 
pfplus_matrix_julia, pfplus_julia = zeros((2^18, n+1)), zeros(n+1)
t1 = time(); for i in 0:(2^18-1) # iterate over 2^m possibilities 
    v = digits(i,base=2,pad=m) # create vector of 0's and 1's from binary number i with in total m digits 
    myset = findall(v.==1) # find places of the 1-elements
    if length(myset)==0
        pfplus_julia .+= ((-1)^length(myset)) .* prod(1e-50 .+ prevminneg .+ (1 .- prev), dims=2)
        pfplus_matrix_julia[i+1,:] = ((-1)^length(myset)) .* prod(1e-50 .+ prevminneg .+ (1 .- prev), dims=2)
    elseif length(myset)==1
        pfplus_julia .+= ((-1)^length(myset)) .* prod(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev), dims=2)
        pfplus_matrix_julia[i+1,:] = ((-1)^length(myset)) .* prod(1e-50 .+ pfmin[myset, :] .* prevminneg .+ (1 .- prev), dims=2)
        t2 = time(); print("$(myset[1]-1): $v --> $(t2-t1) \n")
    else
        pfplus_julia .+= ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
        pfplus_matrix_julia[i+1,:] = ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
    end     
end; 
pfpluspd = pfplus_julia.*prevend
posterior = pfpluspd[2:end] / pfpluspd[1] 
pfplus2 = sum(pfplus_matrix_julia,dims=1)'
pfpluspd2 = pfplus2.*prevend
posterior2 = pfpluspd2[2:end] / pfpluspd2[1] 

# myset = [1,2,4,7,10]
myset = [1,2]
x_julia = ((-1)^length(myset)) .* prod(1e-50 .+ (prod1(pfmin[myset, :],dims=1) .* prevminneg .+ (1 .- prev)),dims=2)
x_matlab = [0.06549453138050265188 0.01708775786397482080 0.00188666328496350376 0.00315957715996739372 0.00467637131842015170 0.00423505357699416676 0.00250060349097871179 0.00418097702716229043 0.00211482226912254470 0.00699837623089968885 0.01658011226869566049 0.00907247888497892362 0.01495998681224277281 0.00919754992716201487 0.00339144181797071441 0.00211466612653891549 0.00245014388516161859 0.00530963673128143770 0.00213102079997943379 0.01196759934680407915 0.00380132453308211537 0.00259633646467418745 0.00140085953047794235 0.00197976738378084121 0.01309914366843020414 0.00507828991744610665 0.00598851582447933965 0.00265058586990244999 0.00177090164552924081 0.00493556236310895100 0.00243251762099873804 0.00238633677961369800 0.01392475577999217375 0.00738165336893128500 0.01912016546400580053 0.00542662738880661348 0.00965248297976941706 0.01045527816323779668 0.01000968873345345518 0.00392862213103524309 0.00205538974348107046 0.00780537827667419828 0.01818673448893423292 0.00182668222768358449 0.00190845666441797181 0.00197909181799309046 0.00625764893086562957 0.00414262264281432466 0.00258050349791668406 0.01201169453770916523 0.00494069991293435810 0.00242984143615578669 0.00524589444131614169 0.00183628495594975673 0.00429670763387092914 0.00539591157140758334 0.00398270215694637286 0.00187559464522636297 0.00347487263403572708 0.00456807189951203424 0.02074902051185462809 0.00403735410752388210 0.00746418526831803611 0.01062614267413054814 0.00670387847420008450 0.01065551696593388178 0.00216256524963866345 0.00299816855078844377 0.01037043638039324318 0.00756432394506582374 0.00586703773507149442 0.01015027568579617012 0.00935038217085253434 0.00908394189749609643 0.00549678906616968450 0.00756158655155163284 0.00436575094615479411 0.01323771089102567661 0.00409876003972844093 0.01225297117539942728 0.00315926052157383327 0.00352357864458155076 0.00266009943424533671 0.01849976260208716083 0.00365804963941016286 0.02208737875914286442 0.01049972869263069569 0.01171233483198696849 0.01234411779874078526 0.00326793218384231927 0.00464219576380456734 0.00359205668036817038 0.00538508876414389714 0.01590865055211706722 0.00769516731542701712 0.01193607253753830623 0.00414568857126547223 0.00520384986245851344 0.00101387717046015863 0.00399381286577738404 0.01378263591508129514]
(x_julia .- x_matlab')
(pfplus_matrix_matlab .- pfplus_matrix_julia)
pfplus_matrix_matlab

pfplus_matlab'
sum(pfplus_matrix_matlab,dims=1)

