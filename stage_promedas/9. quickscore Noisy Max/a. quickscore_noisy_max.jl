function loop_term_fast(sign::Int64,pfminneg::Vector{T},previn_pfminneg::Matrix{T},one_min_previn::Matrix{T},prod_pfmin::Matrix{T}) where {T <: Real}
    A = prod_pfmin .* previn_pfminneg .+ one_min_previn; 
    B = prod_pfmin .* pfminneg' 
    term = sign .* prod(A) .* [1 ; B' ./ A']; # --> prod(A + (B-A) * I0), with 'I0' the diagonal.  
    return term
end

function quickscore_noisy_max(patientprev::Matrix{Float64}, sens_normal::Matrix{Float64}, sens_medium::Array{Float64,3}, sens_extreme::Matrix{Float64}, method::String)
    float_df = DataFrame("Fl32"=>Float32,"Fl128"=>Float128,"MF3"=>Float64x3,"MF4"=>Float64x4,"BF"=>BigFloat)
    float_type = float_df[1,names(float_df)[occursin.(names(float_df),method)][1]]
    patientprev,sens_normal,sens_medium,sens_extreme = [float_type.(x)  for x in [patientprev,sens_normal,sens_medium,sens_extreme]]; 

    previn, pfmin, pfmin1, pfmin2, pfminneg = patientprev[:,1], sens_extreme', sens_medium[:,:,1]', sens_medium[:,:,2]', sens_normal[:,1]
    (m,n),m1 = (!isempty(pfmin) ? size(pfmin) : (0, length(previn))), size(pfmin1,1) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    println("m = $(m+m1) (postive tests)")
    myset_matrix, myset_matrix1 = collect(powerset([1:m...])), collect(powerset([1:m1...]))
    power_m1 = 2^m1

    threading = occursin("thread",method) ? true : false
    pfplus_matrix = threading ? [zeros(float_type,n+1,1) for i=1:Threads.nthreads()] : zeros(float_type,2^m,2^m1,n+1,1) 
    previn_pfminneg = previn' .* pfminneg' 
    one_min_previn = Array(patientprev[:,2]')
    # one_min_previn = 1 .- previn'

    println("Method: '$method'")
    dt = @elapsed begin
        if threading 
            println("Multi-threading: $(Threads.nthreads()) threads");
            Threads.@threads for idx in ProgressBar(CartesianIndices((1:2^m, 1:2^m1))) # In VS Code --> settings (left below) --> choose settings --> type in "Julia: Num Threads" and click "Edit in settings.json" --> set '"julia.NumThreads": 4' (if you want to use 4 threads, which is maximum for my macbook)
                prod_pfmin = prod(pfmin[myset_matrix[idx[1]],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-idx[2]+1],:],dims=1).*prod(pfmin2[myset_matrix1[idx[2]],:],dims=1)
                sign = (-1)^(length(myset_matrix[idx[1]])+length(myset_matrix1[power_m1-idx[2]+1]))
                pfplus_matrix[Threads.threadid()] = pfplus_matrix[Threads.threadid()] .+ loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
            end 
            pfplus = sum(cat(pfplus_matrix...,dims=3),dims=3)
        else
            for idx in ProgressBar(CartesianIndices((1:2^m, 1:2^m1))) # iterate over 2^m possibilities 
                prod_pfmin = prod(pfmin[myset_matrix[idx[1]],:],dims=1) .* prod(pfmin1[myset_matrix1[power_m1-idx[2]+1],:],dims=1).*prod(pfmin2[myset_matrix1[idx[2]],:],dims=1)
                sign = (-1)^(length(myset_matrix[idx[1]])+length(myset_matrix1[power_m1-idx[2]+1]))
                pfplus_matrix[idx[1],idx[2],:,:] = loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin)
                # println(vcat(digits(Int(sum(2 .^myset_matrix1[idx[2]]./2)),base=2,pad=4),digits(Int(sum(2 .^myset_matrix[idx[1]]./2)),base=2,pad=3)))
                # println(Float64.((prod(pfmin1[myset_matrix1[power_m1-idx[2]+1],:],dims=1).*prod(pfmin2[myset_matrix1[idx[2]],:],dims=1))[1:5]'))
                # println(Float64.((loop_term_fast(sign,pfminneg,previn_pfminneg,one_min_previn,prod_pfmin).*[1;previn])[1:6]),'\n')
            end    
            pfplus = sum(pfplus_matrix,dims=(1,2))[1,1,:,:]
        end
    end
    P_joint = pfplus[2:end,1] .* previn
    posterior = P_joint / pfplus[1,1] 
    println("Running time: $dt\n")
    return pfplus_matrix, pfplus, P_joint, posterior, dt 
end 


matlab_dir = "/Users/sam/Documents/MATLAB/Wim/variables/noisy_max/patient_input.mat"; 
matfile = matopen(matlab_dir,"r"); 
    patientprev = read(matfile, "patientprev")[:,1]; 
    sens_normal = read(matfile, "sens_normal")[:,1]; 
    sens_medium = read(matfile, "sens_medium"); 
    sens_extreme = read(matfile, "sens_extreme"); 
close(matfile)


# previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(13,n_disease=1000);
# pfplus_matrix, pfplus, P_joint, posterior, dt = quickscore_noisy_max(previn, pfmin, pfminneg,"trick Fl128 thread");  
# pfplus_matrix_, pfplus_, P_joint_, posterior_, dt_ = quickscore(previn, pfmin, pfminneg,"trick Fl128 thread");  



