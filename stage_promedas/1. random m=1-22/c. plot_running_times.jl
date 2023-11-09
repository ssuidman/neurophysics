include("../packages.jl")
include("../run_one_time.jl")
matlab_dir = "/Users/sam/Documents/MATLAB/Wim/"
m_variables_dir = "variables/random_m_1_22/"
figure_dir = "figures/"



function plotting_m_1_22(save_figure,matlab_dir,m_variables_dir,figure_dir)
    """
    Plotting the running times for Julia and MATLAB
    """
    #     dt_julia, dt_matlab, posterior_julia, posterior_matlab = Vector{Float64}(), Vector{Float64}(), Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    dt_julia = zeros(22)
    dt_matlab = zeros(22)
    for i in 1:22
        previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=i,n_myset=i,language="Julia");
        dt_julia[i] = dt
        previn, pfmin, pfminneg, pfplus, P_joint, posterior, dt, prev, prevminneg, myset = run_one_time_var(m=i,n_myset=i,language="MATLAB");
        dt_matlab[i] = dt
    end

    gr()
    plot!()
    p1 = plot(1:22,dt_julia,label="Julia",legend=:topleft)
    plot!(1:22,dt_matlab,label="MATLAB")
    ylabel!("t(s)")

    p2 = plot(1:22,dt_julia,label="Julia",yscale= :log10,legend=:topleft)
    plot!(1:22,dt_matlab,label="MATLAB",yscale= :log10)
    xlabel!("postive tests 'm'")

    p = plot(p1,p2,layout=(1,2),suptitle="Running time quickscore",dpi=1200)
    xlabel!("postive tests 'm'")
    if save_figure
        savefig(figure_dir*"running_time_comparison_quickscore.png")
    end
end
plotting_m_1_22(true,matlab_dir,m_variables_dir,figure_dir)


