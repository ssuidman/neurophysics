include("packages.jl")
include("quickscore_algorithm.jl")

function run_one_time_var(;m=17,n_myset=0,language="Julia")
    """
    Function with which you can run the quickscore algorithm one time easily and also the important line of the function 'loop_term()'
    You can define m and the length of myset(=0,1,..,m), but this is not necessary. 
    """
    Random.seed!(1234)
    # LOAD SOME RANDOM VARIABLES previn,pfminneg,pfmin FOR A CERTAIN m=17 and run the quickscore
    if language=="Julia"
        file = jldopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/variables/random_m_1_22/m_$m.jld", "r")
            previn = file["previn"]
            pfminneg = file["pfminneg"]
            pfmin = file["pfmin"]
            pfplus = file["pfplus"]
            P_joint = file["P_joint"]
            posterior = file["posterior"]
            dt = file["dt"]
        close(file)
    elseif language=="MATLAB"
        matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
        matfile = matopen(joinpath(matlab_dir,"variables//random_m_1_22/m_$m.mat"),"r")
            previn = read(matfile, "previn")
            pfminneg = read(matfile, "pfminneg")
            pfmin = read(matfile, "pfmin")
            pfplus = read(matfile, "pfplus")
            P_joint = read(matfile, "P_joint")[1,:]
            posterior = read(matfile, "posterior")[1,:]
            dt = read(matfile, "dt")
        close(matfile)
    end

    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    myset = sort(sample(1:m,n_myset,replace=false))

    return previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset
end
# # Let these lines be outcommented so they are not run when adding 'include("run_one_time.jl")' to the header of a file
previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");
# pfplus_, P_joint_, posterior_,dt_ = quickscore(previn, pfmin, pfminneg,"prod");  
