include("../0. useful/packages.jl")

function test_giel(method::String)
    """
    This is a stripped down version of the code of Bart Kappen that should have shown that vectorization should show that 
    loops are faster in julia than vectorization. However, when assigning z as Float64 and not Int64 this is not the case. 
    It is also possible to use another package, such as 'Octavian', that is faster in any case. 
    can be faster than the loops. 
    """
    Random.seed!(1234)
    M=20000
    n=100
    x=rand([-1 1],M,n)      # inputs x(mu,i)=-1,1 mu=1:M, i=1:n
    y=rand([-1 1],M,1)        # output y(mu)=-1,1
    # BECAUSE OF THE LINE BELOW THE CODE FOR VECTORIZATION IS SLOW!!! IT SHOULD BE Float64 AND NOT Int64!!!!!
    z=x.*(y*ones(Int, 1,n))
    # z=x.*(y*ones(1,n))
    iter_max=1e3
    m=zeros(1,n)
    c = Matrix{Float64}(undef,M,1)
    
    t = @elapsed begin
        if method=="loop"
            println("loop")
            for k=1:iter_max
                c=zeros(1,M)
                for i=1:M
                    for j=1:n
                        c[i]+=z[i,j]*m[j]
                    end
                end
            end
        elseif method=="loop @inbounds"
            println("loop @inbounds")
            # This is used in the original code of Bert and is in some cases faster than loops 
            for k=1:iter_max
                c=zeros(1,M)
                for i=1:M
                    for j=1:n
                        @inbounds c[i]+=z[i,j]*m[j]
                    end
                end
            end
        elseif method=="vectorization"
            println("vectorization (slow because 'typeof(z)' is Vector{Int64} and not Vector{Float64})")
            for k=1:iter_max
                c=z*m'
            end
        elseif method=="vectorization Octavian"
            println("vectorization Octavian")
            # This does the same as c=z*m', but then with another package "Octavian" that is 4 times faster 
            for k=1:iter_max
                Octavian.matmul!(c, z, m') 
            end
        end
    end
    println("Running time: $t")
end 

# Run this function a few times after eachother to see the difference 
methods = ["loop","loop @inbounds","vectorization","vectorization Octavian"]
for method in methods
    for i=1:5
        test_giel(method)
    end
end




