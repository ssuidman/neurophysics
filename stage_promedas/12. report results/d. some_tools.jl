include("../0. useful/packages.jl")
include("0. quickscore_algorithm_report.jl")

# function to prepare data for patient 404 case with m=9
function prepare_patient_data(case)
    patient_cases_raw = Dict{String,Matrix{Any}}()
    case_nr = parse(Int,case[end])
    file = XLSX.openxlsx("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/patient_cases.xlsx")
        for case_nr=1:3
            patient_cases_raw["case $case_nr"] = file["Sheet$case_nr"][1:end, 1:end]
        end
    close(file)
    data_alisa = DataFrame(Float64.(patient_cases_raw[case][8:end,:]),patient_cases_raw[case][6,:])
    idx_sens = names(data_alisa)[startswith.(names(data_alisa),"Sensitivity")]
    all_sens = Matrix(data_alisa[:,idx_sens])
    pos_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'true'"
    neg_idx = patient_cases_raw[case][4,1:length(idx_sens)] .== "'false'"
    sens = all_sens[:,pos_idx]
    sensneg = all_sens[:,neg_idx]
    previn = data_alisa[:,"prevalence"]
    pfmin = 1 .- sens'
    pfminneg = prod(1 .- sensneg', dims=1)[1,:]

    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    println("n = $n \t m+ = $m \t m- = $(size(sensneg,2))")
    return previn, pfmin, pfminneg
end 
# prod vs exp-sum-log analysis
previn, pfmin, pfminneg = prepare_patient_data("case 3"); 
m,n = size(pfmin)
dt_prod_quickscore = zeros(9) 
dt_exp_sum_log_quickscore = zeros(9) 
error_prod = zeros(9)
error_exp_sum_log = zeros(9)
for i=1:9
    pfplus_matrix, pfplus, P_joint, posterior_prod, dt, dt_array_not = quickscore_report(previn,pfmin[1:i,:],pfminneg,"prod")
    dt_prod_quickscore[i] = dt
    pfplus_matrix, pfplus, P_joint, posterior_exp_sum_log, dt, dt_array_not = quickscore_report(previn,pfmin[1:i,:],pfminneg,"exp-sum-log")
    dt_exp_sum_log_quickscore[i] = dt
    pfplus_matrix, pfplus, P_joint, posterior_BF, dt, dt_array_not = quickscore_report(previn,pfmin[1:i,:],pfminneg,"trick BF thread")
    error_prod[i] = maximum(abs,(posterior_prod .- posterior_BF)./posterior_BF)
    error_exp_sum_log[i] = maximum(abs,(posterior_exp_sum_log .- posterior_BF)./posterior_BF)
end
# writedlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/some_tools.csv",hcat(dt_prod_quickscore,dt_exp_sum_log_quickscore,error_prod,error_exp_sum_log),',')
dt_matrix = DataFrame(readdlm("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/variables/some_tools.csv",','),["prod QS","exp-sum-log QS","prod error","exp-sum-log error"])
# prod vs exp-sum-log figure 
gr() 
default(fontfamily="Times New Roman") 
p1 = plot(1:9,dt_matrix[:,"prod QS"],ylims = (10^-3,10^0.5),label="prod",legend=:topleft,yscale=:log10) 
# p1 = plot(1:9,dt_matrix[:,"prod QS"],ylims = (10^-3,10^0.5),label="prod",legend=:topleft,yscale=:log10,title="running time") 
plot!(1:9,dt_matrix[:,"exp-sum-log QS"],label="exp-sum-log",yscale=:log10) 
ylabel!("t (s)") 
p2 = plot(1:9,dt_matrix[:,"prod error"],label="prod",legend=:topleft,yscale=:log10) 
# p2 = plot(1:9,dt_matrix[:,"prod error"],label="prod",legend=:topleft,yscale=:log10,title="error") 
plot!(1:9,dt_matrix[:,"exp-sum-log error"],label="exp-sum-log",yscale=:log10) 
ylabel!(L"$\delta_{64}$")
p = plot(p1,p2,layout=(1,2),dpi=600) 
xticks!([2,4,6,8], [L"2^2",L"2^{4}",L"2^{6}",L"2^{8}"])
xlabel!(L"#iterations ($2^k$)")
# savefig("/Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12. report results/figures/d_some_tools.png")
# prod vs exp-sum-log profiling tools to show that log is intensive 
@pprof pfplus_matrix, pfplus, P_joint, posterior_prod, dt, dt_array_not = quickscore_report(previn,pfmin,pfminneg,"prod");
PProf.clear()
@pprof pfplus_matrix, pfplus, P_joint, posterior_prod, dt, dt_array_not = quickscore_report(previn,pfmin,pfminneg,"exp-sum-log")
PProf.clear()

# precompilation time 
# Go to the terminal and find the file "d. some_tools_precompilation_script.jl" Then run the file using
# 'julia --trace-compile=stderr d.\ some_tools_precompilation_script.jl'
# Or you can run the line(s) below:
run(`/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --trace-compile=stderr /Users/sam/Documents/Programmeren/neurophysics/stage_promedas/12.\ report\ results/d.\ some_tools_precompilation_script.jl`)
# run(`/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --trace-compile=stderr 12.\ report\ results/d.\ some_tools_precompilation_script.jl`)


