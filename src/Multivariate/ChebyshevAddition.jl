#
# A new ProductFun constructor for bivariate functions defined as the
# difference of their arguments.
#
# This method takes as input a 1D Fun which is the antidiagonal of the
# bivariate function, then forms the matrix of ProductFun coefficients
# recursively via a Chebyshev addition theorem.
#
function ProductFun{S,T}(f::Fun{S,T})
    ## TODO: generalize for more intervals & spaces.
    @assert S == typeof(Ultraspherical{0}(Interval(T[-2,2])))

    c = chop(coefficients(f),10eps(T))
    N = length(c)
    ## TODO: Can the coefficients be stored in only three 2D arrays C_old, C_new, and temp?
    ## TODO: Can the complexity be reduced from O(n^3), albeit with a low constant?
    C,X = zeros(T,N,N,N),zeros(T,N,N)
    un = one(T)

    C[1,1,1] = un
    C[2,1,2] = -un/2
    C[1,2,2] = un/2
    C[1,1,3] = -un/2
    C[3,1,3] = un/4
    C[2,2,3] = -un
    C[1,3,3] = un/4
    for n=4:N
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
        C[1,1,n] = (C[1,2,n-1]-C[2,1,n-1])/2 - C[1,1,n-2]
        C[2,1,n] = (C[2,2,n-1]-C[3,1,n-1])/2 - C[1,1,n-1] - C[2,1,n-2]
        C[n,1,n] = C[n-1,1,n-1]/(-2)
        C[1,2,n] = (C[1,3,n-1]-C[2,2,n-1])/2 + C[1,1,n-1] - C[1,2,n-2]
        C[2,2,n] = (C[2,3,n-1]-C[3,2,n-1])/2 + C[2,1,n-1]-C[1,2,n-1]-C[2,2,n-2]
        C[1,n,n] = C[1,n-1,n-1]/2
        for i=3:n-1
            C[i,1,n] = (C[i,2,n-1]-C[i-1,1,n-1]-C[i+1,1,n-1])/2 - C[i,1,n-2]
            C[i,2,n] = (C[i,3,n-1]-C[i-1,2,n-1]-C[i+1,2,n-1])/2 + C[i,1,n-1] - C[i,2,n-2]
        end
        for j=3:n-1
            C[1,j,n] = (C[1,j+1,n-1]+C[1,j-1,n-1]-C[2,j,n-1])/2 - C[1,j,n-2]
            C[2,j,n] = (C[2,j+1,n-1]+C[2,j-1,n-1]-C[3,j,n-1])/2 - C[1,j,n-1] - C[2,j,n-2]
        end
        for i=3:n
            for j=3:n-i+1
                C[i,j,n] = (C[i,j+1,n-1]+C[i,j-1,n-1]-C[i+1,j,n-1]-C[i-1,j,n-1])/2 - C[i,j,n-2]
            end
        end
    end
    for n=1:N
        for i=1:n
            for j=1:n-i+1
                X[i,j] += c[n]*C[i,j,n]
            end
        end
    end
    V = Ultraspherical{0}()
    ProductFun(X,V⊗V)
end
