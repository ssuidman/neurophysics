include("../0. useful/packages.jl")

function quickscore_preparation(m_wanted;n_disease=100,n_test=30,expected_number_of_diseases=3,sparseness=0.5)
    """
    Function that creates some variables that can be put into the quickscore algorithm
    """
    # the input parameters are 'typical number of diags', 'number of findings', 'typical number of diags', 'connection density'
    # the model
    previn = (rand(n_disease) .+ 1) / 2  # prevalences of diseases, random between 0.5 and 1
    previn .= expected_number_of_diseases * previn / sum(previn)  # in expectation sumofdiags diseases present

    connections = rand(n_test, n_disease) .< sparseness  # selection of connections
    sensorig = (1 .+ 4 .* rand(n_test, n_disease)) ./ 10 .* connections  # sensitivities random between 0.1 and 0.5

    # data generation. Generate diseases according to the prevalences. Then test results (=findings) according to the sensitivities.
    # accept if the number of diseases and the number of positive findings is OK, else try again.
    actualdiags, postest, negtest, ndiags,m,testresults = 0,0,0,0,0,0 #setting some variables to use outside the loop too 

    # while (ndiags==0) | !(m>5 && m<15) # changing the m range, could change the amount postest/negtest which both should never be empty 
    while (ndiags==0) | !(m==m_wanted) # changing the m range, could change the amount postest/negtest which both should never be empty 
        actualdiags = findall(rand(n_disease) .< previn)  # generate diseases and get the y value out of a 1D-array of cartesian coordinates
        ndiags = length(actualdiags) 

        testresults = vec( any(rand(n_test, ndiags).<(sensorig[:,actualdiags]), dims=2) )  # generate test results

        postest = findall(testresults.==1)
        negtest = findall(testresults.==0)
        m = length(postest)
    end

    actualdiseases = zeros(size(previn))
    actualdiseases[actualdiags] .= 1

    # gather data (some data preparation for the Quickscore algorithm)
    sens = sensorig[postest, :]
    sensneg = sensorig[negtest, :]

    pfmin = 1 .- sens 
    pfminneg = vec( prod(1 .- sensneg, dims=1) )  # be careful with the dimensions: scalar result prod(1 - sensneg) for singleton neg test is wrong.

    m,n = !isempty(pfmin) ? size(pfmin) : (0, length(previn)) # m number of postests. n number of diags. if size(postest) is 0, then we set m=0 `` 
    # println("m = $m (postive tests)")
    prev = repeat(previn',inner=(n+1,1)) # copies of prevalences
    for i in 1:n prev[i+1,i]=1 end # set the (i+1,i) entry 1 to condition on d_i=1,  etc. 
    prevminneg = !isempty(pfminneg) ? prev.*pfminneg' : prev.*ones(n+1,n) # absorb negative findings in prevminneg (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) ), which is needed in heckerman eqn 11 (compare with heckerman 10, there is a constant difference which is incorporated in the block below)
    
    return previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg
end 


previn, pfmin, pfminneg, sens, sensneg, prev, prevminneg = quickscore_preparation(20;n_disease=1000,n_test=30,expected_number_of_diseases=3,sparseness=0.5)
quickscore(previn,pfmin,pfminneg,"MATLAB")
