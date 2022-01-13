
using .LoopVectorization: @avx
function InnerProductV(A::AbstractMatrix, x::AbstractVector)
    @boundscheck @assert size(A,1) == size(A,2) == length(x)
    s = zero(suff(x))
    @inbounds @avx for j in eachindex(x)
        for i in eachindex(x)
            s += x[j] * A[i,j] * x[i]
        end
    end;    s
end
function InnerProductV(A::Diagonal, x::AbstractVector)
    @boundscheck @assert size(A,1) == size(A,2) == length(x)
    s = zero(suff(x))
    @inbounds @avx for j in eachindex(x)
        s += A.diag[j] * x[j] * x[j]
    end;    s
end

@info "InformationGeometry: Loaded LoopVectorization extensions."
