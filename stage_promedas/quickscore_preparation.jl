using Random

# parameters
n_disease = 40   # number of diseases
n_test = 30      # number of findings
expected_number_of_diseases = 3  # typical number of diags
sparseness1 = 0.5  # connection density

# the model
prev = (rand(1, n_disease) .+ 1) / 2  # prevalences of diseases, random between 0.5 and 1
prev .= expected_number_of_diseases * prev / sum(prev)  # in expectation sumofdiags diseases present


connections = rand(n_test, n_disease) .< sparseness1  # selection of connections
sensorig = (1 .+ 4 .* rand(n_test, n_disease)) ./ 10 .* connections  # sensitivities random between 0.1 and 0.5

#=
# data generation. Generate diseases according to the prevalences. Then
# test results (=findings) according to the sensitivities.
# accept if the number of diseases and the number of positive findings is OK,
# else try again.
flag = true
while flag
    actualdiags = find(rand(1, n_disease) .< prev)  # generate diseases
    ndiags = length(actualdiags)

    testresults = any(rand(n_test, ndiags) .< (sensorig[:, actualdiags]), dims=2)  # generate test results

    postest = findall(testresults .== 1)
    negtest = findall(testresults .== 0)
    m = length(postest)
    println([ndiags, m])
    if (m > 5) && (m < 15) && (!isempty(actualdiags))
        # for now, only cases with at least 5 positive tests, but no more than
        # 14 positive tests, and at least one disease.
        flag = false
    end
end

# gather data (some data preparation for the Quickscore algorithm)
sens = sensorig[postest, :]
sensneg = sensorig[negtest, :]

pfmin = 1 .- sens
pfminneg = prod(1 .- sensneg, dims=1)  # be careful with the dimensions: scalar result prod(1 - sensneg) for singleton neg test is wrong.

# Here it goes into quickscore.
tic()
(posterior, pfpluspd, pfplus) = fun_quickscoreW(prev, pfmin, pfminneg)
timep = toq()

actualdiseases = zeros(size(prev))
actualdiseases[actualdiags] .= 1
println("prevalences, posteriors, actual diseases")
println([prev', posterior', actualdiseases'])
println("\nPositive findings: $m, time: $timep")
```

Please note that there may be some differences in indexing and broadcasting between MATLAB and Julia, so I've adjusted the code accordingly. Also, Julia uses 1-based indexing, so some array indices have been adjusted accordingly.
=#
