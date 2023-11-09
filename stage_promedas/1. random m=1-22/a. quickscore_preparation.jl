include("../packages.jl")

function quickscore_preparation(m_wanted)
    """
    Function that creates some variables that can be put into the quickscore algorithm
    """
    # the input parameters are 'typical number of diags', 'number of findings', 'typical number of diags', 'connection density'
    n_disease, n_test, expected_number_of_diseases, sparseness = 100, 30, 3, 0.5;
    # the model
    prev = (rand(n_disease) .+ 1) / 2  # prevalences of diseases, random between 0.5 and 1
    prev .= expected_number_of_diseases * prev / sum(prev)  # in expectation sumofdiags diseases present

    connections = rand(n_test, n_disease) .< sparseness  # selection of connections
    sensorig = (1 .+ 4 .* rand(n_test, n_disease)) ./ 10 .* connections  # sensitivities random between 0.1 and 0.5

    # data generation. Generate diseases according to the prevalences. Then test results (=findings) according to the sensitivities.
    # accept if the number of diseases and the number of positive findings is OK, else try again.
    actualdiags, postest, negtest, ndiags,m,testresults = 0,0,0,0,0,0 #setting some variables to use outside the loop too 

    # while (ndiags==0) | !(m>5 && m<15) # changing the m range, could change the amount postest/negtest which both should never be empty 
    while (ndiags==0) | !(m==m_wanted) # changing the m range, could change the amount postest/negtest which both should never be empty 
        actualdiags = findall(rand(n_disease) .< prev)  # generate diseases and get the y value out of a 1D-array of cartesian coordinates
        ndiags = length(actualdiags) 

        testresults = vec( any(rand(n_test, ndiags).<(sensorig[:,actualdiags]), dims=2) )  # generate test results

        postest = findall(testresults.==1)
        negtest = findall(testresults.==0)
        m = length(postest)
    end

    actualdiseases = zeros(size(prev))
    actualdiseases[actualdiags] .= 1

    # gather data (some data preparation for the Quickscore algorithm)
    sens = sensorig[postest, :]
    sensneg = sensorig[negtest, :]

    pfmin = 1 .- sens 
    pfminneg = vec( prod(1 .- sensneg, dims=1) )  # be careful with the dimensions: scalar result prod(1 - sensneg) for singleton neg test is wrong.

    return prev, pfmin, pfminneg
end 



