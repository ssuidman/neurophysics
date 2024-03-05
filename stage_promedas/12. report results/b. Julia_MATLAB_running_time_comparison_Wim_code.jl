include("../0. useful/packages.jl")

figure_dir = "figures/"
function plotting_m_1_22(save_figure)
    """
    Plotting the running times for Julia and MATLAB
    """
    dt_julia = zeros(22)
    dt_matlab = zeros(22)
    for i in 1:22
        file = jldopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/m_1_22_julia/m_$i.jld", "r")
            pfplus_julia = file["pfplus"]
            posterior_julia = file["posterior"]
            dt_julia[i] = file["dt"]
        close(file)
        matfile = matopen("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/m_1_22_matlab/m_$i.mat","r")
            pfplus_matlab = read(matfile, "pfplus")[1,:]
            posterior_matlab = read(matfile, "posterior")[1,:]
            dt_matlab[i] = read(matfile, "dt")
        close(matfile)
    end

    gr()
    plot!()
    default(fontfamily="Times New Roman")
    p1 = plot(1:22,dt_julia,label="Julia",legend=:topleft)
    plot!(1:22,dt_matlab,label="MATLAB")
    
    p2 = plot(1:22,dt_julia,label="Julia",ylims=(10^-2.99,10^2.99),yscale= :log10,legend=:topleft)
    plot!(1:22,dt_matlab,label="MATLAB",yscale= :log10)
    
    p = plot(p1,p2,layout=(1,2),dpi=600)
    # p = plot(p1,p2,layout=(1,2),suptitle="Running time quickscore",dpi=1200)
    # xticks!([5,10,15,20], [L"2^5",L"2^{10}",L"2^{15}",L"2^{20}"])
    # xticks!([5,10,15,20], [L"2^5",L"2^{10}",L"2^{15}",L"2^{20}"])
    # xlabel!(L"#iterations ($2^{m^+}$)")
    xticks!([5,10,15,20], [L"2^{5}",L"2^{10}",L"2^{15}",L"2^{20}"])
    xlabel!(L"#iterations ($2^k$)")
    ylabel!("t (s)")
    
    display(p)
    if save_figure
        savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/b_running_time_comparison_quickscore.png")
    end
end
# plotting_m_1_22(true)


