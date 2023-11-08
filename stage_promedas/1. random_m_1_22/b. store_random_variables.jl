include("../packages.jl")
include("a. quickscore_preparation.jl")
include("../quickscore_algorithm.jl")
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
m_variables_dir = "variables/random_m_1_22/"



function store_random_pfmin(matlab_dir,m_variables_dir)
    """ 
    This function creates for m=1-22 random variables for previn, pfmin, pfminneg
    It then stores it as julia variables (JDL-file) and matlab variables (MAT-file)
    IT ALSO REMOVES 'posterior', etc SO BE CAREFUL!!!
    """
    Random.seed!(1234)
    for i in 1:22
        previn, pfmin, pfminneg = quickscore_preparation(i)

        # Write variables to jld-file 
        jldopen(m_variables_dir*"m_$i.jld","w") do file
            file["previn"] = previn
            file["pfmin"] = pfmin
            file["pfminneg"] = pfminneg
        end
        
        # Write variables to mat-file
        matfile = matopen(matlab_dir*m_variables_dir*"m_$i.mat","w")
            write(matfile,"previn",previn)
            write(matfile,"pfmin",pfmin)
            write(matfile,"pfminneg",pfminneg)
        close(matfile)
    end
    println("Creating 'previn', 'pfmin', 'pfminneg'")
end
################## DO NOT RUN THIS, BECAUSE ALL VARIABLES WILL RESTORE!!! ##################
# store_random_pfmin(matlab_dir,m_variables_dir)



function store_random_posterior(m_variables_dir,run_julia,run_matlab)
    """
    This function runs the quickscore algorithm for all m=1-22 in Julia and MATLAB
    It then stores in the variables in the same JDL- and MAT-files as before 
    """
    # RUNNING THE QUICKSCORE ALGORITHM FROM JL-FILES FOR ALL M=1,...,22
    if run_julia
        for i in 1:22
            jldopen(m_variables_dir*"m_$i.jld", "r") do file
                global previn = file["previn"]
                global pfminneg = file["pfminneg"] 
                global pfmin = file["pfmin"] 
            end
            
            pfplus, P_joint, posterior,dt = quickscore(previn, pfmin, pfminneg,"MATLAB");  

            jldopen(m_variables_dir*"m_$i.jld","a") do file # make sure that this part has not yet run, otherwise it gives an error that parts of the file already exist
                file["pfplus"] = pfplus 
                file["P_joint"] = P_joint 
                file["posterior"] = posterior
                file["dt"] = dt
            end
        end
    end
    if run_matlab
        run(`/Applications/MATLAB_R2022b.app/bin/matlab -nodesktop -nosplash -r "run('~/Documents/MATLAB/Wim/run_all.m');exit;"`)
    end
end
################## DO NOT RUN THIS, BECAUSE ALL VARIABLES WILL RESTORE!!! ##################
# store_random_posterior(m_variables_dir,true,true)





