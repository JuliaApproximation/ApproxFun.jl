


function defaultaxpy!(a,X::SubBandedOperator,Y::BandedMatrix)
     @assert size(X)==size(Y)
     @assert bandwidth(X,1) ≤ bandwidth(Y,1) && bandwidth(X,2) ≤ bandwidth(Y,2)

     br=bandinds(X)
     for (k,j) in eachbandedindex(X)
         Y[k,j]+=a*X[k,j]
     end

     Y
end

function defaultaxpy!(a,X::SubBandedOperator,Y)
     @assert size(X)==size(Y)

     br=bandinds(X)
     for (k,j) in eachbandedindex(X)
         Y[k,j]+=a*X[k,j]
     end

     Y
end

BLAS.axpy!(a,X::SubBandedOperator,Y::AbstractMatrix) = defaultaxpy!(a,X,Y)

