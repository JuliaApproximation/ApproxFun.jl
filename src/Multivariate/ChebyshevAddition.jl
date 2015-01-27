#
# A new ProductFun constructor for bivariate functions defined as the
# difference of their arguments.
#
# This method takes as input a 1D Fun which is the antidiagonal of the
# bivariate function, then forms the matrix of ProductFun coefficients
# recursively via a Chebyshev addition theorem.
#
function ProductFun{S,T}(f::Fun{S,T})
    ## TODO: Can the complexity be reduced from O(n^3)?
    ## TODO: generalize for more intervals & spaces.
    @assert S == typeof(Ultraspherical{0}(Interval(T[-2,2])))
    c = chop(coefficients(f),10eps(T))
    N = length(c)
    if N ≤ 3 N=3;pad!(c,3) end
    un = one(T)
    C1,C2,X = zeros(T,N,N),zeros(T,N,N),zeros(T,N,N)

    C1[1,1] = un
    cn = c[1]

    X[1,1] += cn*C1[1,1]

    C2[2,1] = -un/2
    C2[1,2] = un/2
    cn = c[2]

    X[2,1] += cn*C2[2,1]
    X[1,2] += cn*C2[1,2]

    C1[1,1] = -un/2
    C1[3,1] = un/4
    C1[2,2] = -un
    C1[1,3] = un/4
    cn = c[3]

    X[1,1] += cn*C1[1,1]
    X[3,1] += cn*C1[3,1]
    X[2,2] += cn*C1[2,2]
    X[1,3] += cn*C1[1,3]

    @inbounds for n=4:N
        #
        # There are 11 unique recurrence relationships for the coefficients. The main recurrence is:
        #
        # C[i,j,n] = (C[i,j+1,n-1]+C[i,j-1,n-1]-C[i+1,j,n-1]-C[i-1,j,n-1])/2 - C[i,j,n-2],
        #
        # and the other 10 come from shutting some terms off if they are out of bounds,
        # or for the row C[2,1:n,n] or column C[1:n,2,n] terms are turned on. This follows from
        # the reflection of Chebyshev polynomials: 2T_m(x)T_n(x) = T_{m+n}(x) + T_|m-n|(x).
        # For testing of stability, they should always be equal to:
        # C[1:n,1:n,n] = coefficients(ProductFun((x,y)->cos((n-1)*acos((y-x)/2)))).
        #
        C2[1,1] = (C1[1,2]-C1[2,1])/2 - C2[1,1]
        C2[2,1] = (C1[2,2]-C1[3,1])/2 - C1[1,1] - C2[2,1]
        C2[n,1] = C1[n-1,1]/(-2)
        C2[1,2] = (C1[1,3]-C1[2,2])/2 + C1[1,1] - C2[1,2]
        C2[2,2] = (C1[2,3]-C1[3,2])/2 + C1[2,1]-C1[1,2] - C2[2,2]
        C2[1,n] = C1[1,n-1]/2
        for k=n-2:-2:3
            C2[k,1] = (C1[k,2]-C1[k-1,1]-C1[k+1,1])/2 - C2[k,1]
            C2[1,k] = (C1[1,k+1]+C1[1,k-1]-C1[2,k])/2 - C2[1,k]
        end
        for k=n-1:-2:3
            C2[k,2] = (C1[k,3]-C1[k-1,2]-C1[k+1,2])/2 + C1[k,1] - C2[k,2]
            C2[2,k] = (C1[2,k+1]+C1[2,k-1]-C1[3,k])/2 - C1[1,k] - C2[2,k]
        end
        for j=n:-1:3,i=n-j+1:-2:3
            C2[i,j] = (C1[i,j+1]+C1[i,j-1]-C1[i+1,j]-C1[i-1,j])/2 - C2[i,j]
        end

        cn = c[n]
        for j=n:-1:1,i=n-j+1:-2:1
            X[i,j] += cn*C2[i,j]
        end

        for j=1:n,i=1:n-j+1
            C1[i,j],C2[i,j] = C2[i,j],C1[i,j]
        end
    end
    V = Ultraspherical{0}()
    ProductFun(X,V⊗V)
end
