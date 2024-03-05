# Import the packages that are needed
# To install a package, you can run: 
# using Pkg
# Pkg.add("MAT")
# Pkg.add("DelimitedFiles")
# etc
using MAT
using DelimitedFiles
using Quadmath
using MultiFloats
using DataFrames
using Combinatorics
using ProgressBars
using LinearAlgebra


function loop_term_noisy_max(sign::Int64,pfminneg::Vector{T},previn_pfminneg::Matrix{T},one_min_previn::Matrix{T},prod_pfmin::Matrix{T}) where {T <: Real}
    """
    Function that does the calculation of the loop terms and returns a (n+1,1)-sized matrix for 1 value of p(F+,F-) and n values of p(F+,F-|di+)
    """
    A = prod_pfmin .* previn_pfminneg .+ one_min_previn; # quickscore line for p(F+,F-)
    B = prod_pfmin .* pfminneg' # quickscore line for p(F+,F-|di+)
    # What you essentially do for p(F+,F-) and all p(F+,F-|di+) is prod(A + (B-A) * I0), with " I0=[zeros(n) I(n)]' " the diagonal.  
    # prod(A) calculates prod for p(F+,F-), but the fraction B'./A' is a very efficient way to calculate the term for all di using a one-time-calculation prod(A)
    term = sign .* prod(A) .* [1 ; B' ./ A']; 
    return term
end

function quickscore_noisy_max(patientprev::Matrix{Float64}, sens_normal::Matrix{Float64}, sens_medium::Array{Float64,3}, sens_extreme::Matrix{Float64}, method::String)
    """
    The entire quickscore algorithm that works for noisy max and takes 4 input variables from MATLAB:
        patientprev     -->     prevalences in a (n,2)-matrix which is [prev,1-prev]
        sens_normal     -->     negative findings
        sens_medium     -->     findings for 'below' and 'above'
        sens_extreme    -->     findings for 'true', 'ebelow', 'eabove'
    from MATLAB. Also a method needs to be chosen that defines the precision (Float32, Float128, Float64x3, Float64x4, BigFloat) and whether or 
    not you want to use threading. 
    The function returns: 
        pfplus          -->     p(F+,F-) and p(F+,F-|di+) put together in a (n+1,1)-sized matrix
        P_joint         -->     p(F+,F-|di+)*p(di+) in a n-sized vector
        posterior       -->     p(di+|F+,F-) which is called 'pdiag' in the MATLAB-code
        dt              -->     running time of the loop-part of the quickscore algorithm 
    """
    # Changing the float type based on whether "Fl32", "Fl128", etc is present in the method, Float64 is default. 
    float_df = DataFrame("Fl32"=>Float32,"Fl128"=>Float128,"MF3"=>Float64x3,"MF4"=>Float64x4,"BF"=>BigFloat)
    float_type = float_df[1,names(float_df)[occursin.(names(float_df),method)][1]]
    patientprev,sens_normal,sens_medium,sens_extreme = [float_type.(x)  for x in [patientprev,sens_normal,sens_medium,sens_extreme]]; 

    # Preparation of quickscore variables based on the MATLAB-input variables 
    previn, pfmin, pfmin1, pfmin2, pfminneg = patientprev[:,1], sens_extreme', sens_medium[:,:,1]', sens_medium[:,:,2]', sens_normal[:,1]
    # Obtaining the needed sizes 
    (m,n),m1 = (!isempty(pfmin) ? size(pfmin) : (0, length(previn))), size(pfmin1,1) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    power_m1 = 2^m1
    println("m = $(m+m1) (postive tests)")
    # Creating the set of all possible combinations needed for quickscore 
    myset_matrix, myset_matrix1 = collect(powerset([1:m...])), collect(powerset([1:m1...]))

    # Setting threading to true or false based on whether "thread" is in the input method
    threading = occursin("thread",method) ? true : false

    # Creating the input matrix (2 possible forms) in which all summing terms are saved in, based on whether threading is chosen or not 
    pfplus_matrix = threading ? [zeros(float_type,n+1,1) for i=1:Threads.nthreads()] : zeros(float_type,2^m,2^m1,n+1,1) 
    # Calculating two other input variables before the loop to spare time 
    previn_pfminneg = previn' .* pfminneg' 
    # one_min_previn = Array(patientprev[:,2]')
    one_min_previn = 1 .- previn'

    println("Method: '$method'")
    dt = @elapsed begin
        if threading 
            println("Multi-threading: $(Threads.nthreads()) threads");
            # If threading it chosen then Threads.@threads is used and all combinations of idx=(i,j) are picked 
            # If you use Visual Studio Code you need to change a file in settings (left below) --> choose settings
            # then you type in "Julia: Num Threads" and click "Edit in settings.json"
            # thereafter you set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
            Threads.@threads for idx in ProgressBar(CartesianIndices((1:2^m, 1:2^m1))) 
            # Threads.@threads for idx in CartesianIndices((1:2^m, 1:2^m1)) 
                # prod_pfmin is the part that changes each iteration due to the changing input set from myset_matrix
                prod_pfmin = prod(pfmin[myset_matrix[idx[1]],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-idx[2]+1],:],dims=1).*prod(pfmin2[myset_matrix1[idx[2]],:],dims=1)
                # calculating the sign 
                sign = (-1)^(length(myset_matrix[idx[1]])+length(myset_matrix1[power_m1-idx[2]+1]))
                # the term that is obtained from the loop-function is stored in pfplus_matrix 
                pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term_noisy_max(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
            end 
            # all terms from pfplus_matrix are summed to get pfplus 
            pfplus = sum(cat(pfplus_matrix...,dims=3),dims=3)
        else
            for idx in ProgressBar(CartesianIndices((1:2^m, 1:2^m1))) 
            # for idx in CartesianIndices((1:2^m, 1:2^m1)) # iterate over 2^m possibilities 
                # prod_pfmin is the part that changes each iteration due to the changing input set from myset_matrix
                prod_pfmin = prod(pfmin[myset_matrix[idx[1]],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-idx[2]+1],:],dims=1).*prod(pfmin2[myset_matrix1[idx[2]],:],dims=1)
                # calculating the sign 
                sign = (-1)^(length(myset_matrix[idx[1]])+length(myset_matrix1[power_m1-idx[2]+1]))
                # the term that is obtained from the loop-function is stored in pfplus_matrix 
                pfplus_matrix[idx[1],idx[2],:,:] = loop_term_noisy_max(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
            end    
            # all terms from pfplus_matrix are summed to get pfplus 
            pfplus = sum(pfplus_matrix,dims=(1,2))[1,1,:,:]
        end
    end
    # calculating P_joint and posterior from pfplus and returing it 
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus, P_joint, posterior, dt 
end 

# Get variables that were stored in the .mat-file to use in the quickscore algorithm 
matfile = matopen("patient_input.mat","r"); 
patientprev = read(matfile, "patientprev"); 
sens_normal = read(matfile, "sens_normal"); 
sens_medium = read(matfile, "sens_medium"); 
sens_extreme = read(matfile, "sens_extreme"); 
close(matfile)
# Get the method stored in the csv-file to use in the quickscore algorithm 
method = String(readdlm("patient_method.csv",',')[1]);

# Run the quickscore algorithm 
pfplus, P_joint, posterior, dt = quickscore_noisy_max(patientprev,sens_normal,sens_medium,sens_extreme,method);

# return pdiag and the other variables to a .mat-file 
matfile = matopen("patient_output.mat","w")
    write(matfile,"pfplus",Float64.(pfplus))
    write(matfile,"P_joint",Float64.(P_joint))
    write(matfile,"pdiag",Float64.(posterior))
    write(matfile,"dt",Float64.(dt))
close(matfile)

