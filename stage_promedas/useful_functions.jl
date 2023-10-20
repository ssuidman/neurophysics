function prod1(x;dims=1)
    if dims==1
        result = zeros(1,size(x,2))
        for i=axes(x,2)
            result[1,i] = prod(x[:,i])
        end
    elseif dims==2
        result = zeros(size(x,1),1)
        for i=axes(x,1)
            result[i,1] = prod(x[i,:])
        end
    end    
    return result
end

struct Patient
    age::Union{Int,Nothing}
    gender::Union{String, Nothing}
    test::Union{Array{Int},Nothing}
    tvalue::Union{Array{Any},Nothing}
    diag::Union{Array{Int},Nothing}
    dvalue::Union{Array{Any},Nothing}
    cdiag::Union{Vector{Int},Nothing}
    ext::Union{Vector{Int},Nothing}
    evalue::Union{Array{Any},Nothing}
end


