# module test

include("packages.jl")

function main(with_loop)
    Random.seed!(1234)
    M=20000
    n=100
    x=rand([-1 1],M,n)      # inputs x(mu,i)=-1,1 mu=1:M, i=1:n
    y=rand([-1 1],M,1)        # output y(mu)=-1,1
    # z=x.*(y*ones(Int, 1,n))
    z=x.*(y*ones(1,n))
    beta=2
    dw_max=1
    iter_max=1e3
    m=zeros(1,n)
    learning=1
    iter=0
    # c = 0
    c = Matrix{Float64}(undef,M,1)
    t1 = time()
    # while (learning==1)&(iter<iter_max)
    for k=1:iter_max
        iter=iter+1
        learning=0
        if with_loop
            c=zeros(1,M)
            for i=1:M
                for j=1:n
                    c[i]+=z[i,j]*m[j]
                    # @inbounds removes array bound checks, only use it if you're sure that you won't go out of bounds
                    # @inbounds c[i]+=z[i,j]*m[j]
                end
            end
        else
            c=z*m'
            # Octavian.matmul!(c, z, m') 
        end
        # for mu=1:M
        #     if c[mu]<=0
        #         m=m+z[mu,:]'
        #         learning=1
        #     end
        # end
    end
    t2 = time()
    scaling=maximum(abs,m)
    m=m./scaling
    sigma_h=sum(1 .- m.^2)

    klad = sigmoid.(sqrt(8/pi) * z.*m /sqrt(sigma_h))
    likelihood=sum(log.(klad))/M
    correct=sum(klad.>0.5)/M
    if with_loop
        println("### with loop ###")
    else
        println("### with  *   ###")
    end
    # println("alpha      = ", M/n)
    # println("iter       = $iter")
    # println("m^2        = ", sum(m.^2))
    # println("likelihood = $likelihood")
    # println("correct    = $correct")
    # println("correct    = %6.4f\n", correct)
    println("Running time: $(t2-t1)")
end # function main

function sigmoid(x) 
    return 1 ./ (1+exp(-2*x)) 
end

# The function and the functions it calls are compiled at first invocation.
# This takes significant time on my machine so don't count the ``warm-up call''.
main(true)
main(true)
main(true)
main(true)
main(true)
main(true)
main(false)
main(false)
main(false)
main(false)
main(false)
main(false)
print("\n\n\n")
# end


