# This file shows some profiling tools that can be helpful
include("../0. useful/packages.jl")
include("../0. useful/quickscore_algorithm.jl")
include("../0. useful/run_one_time.jl")

println("Do not run this file in the terminal!!!")
# Get variables for a random case m=17
previn, pfmin, pfminneg, pfplus_, P_joint_, posterior_, dt_, prev, prevminneg, myset = run_one_time_var(m=17,n_myset=3,language="Julia");

# With @pprof it opens a link to localhost where you can see what core functions are used the most 
println("########## @pprof ############")
@pprof pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"exp-sum-log");
PProf.clear()
println("########## @pprof ############")
@pprof pfplus, P_joint, posterior, dt = quickscore(previn,pfmin,pfminneg,"prod");
PProf.clear()
# Run this to see what can cause errors when you will compile the code --> not possible in the terminal, only jl-file in VS code
@code_warntype quickscore(previn, pfmin, pfminneg, "MATLAB")           
# Run this to get information about how much running time a code line or function will take      
@benchmark log.(rand(101,100))
# Run the line below or in your mac terminal this: "julia --trace-compile=stderr e.\ precompile_test.jl" to see 
# that precompiling does not happend during the loops
run(`/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --trace-compile=stderr "2. running time comparison/d. precompile_test.jl"`)
