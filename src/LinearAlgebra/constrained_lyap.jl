
# permute the columns until the principle subblock is nonsingular
# assume B approximates K x ∞
function nonsingular_permute(B)
    K = size(B,1)

    kr=[1]
    while rank(B[:,kr])<K
        if rank(B[:,kr])==length(kr)
            kr=[kr;kr[end]+1]  # add another col
        else
            kr[end]+=1  # integrate last col by 1
        end
    end

    kr=[kr;setdiff(1:size(B,2),kr)]
    eye(size(B,2))[:,kr]
end

## rearrange the rows and columns of bcs so that the principle block
# is identity
# assume we have the condition
#    BX=G
# where B is  K x n, X is n x m and G is K x *
#
# We also need to update L & M in
#    LX* + MX* = *
# to reflect the new ordering of rows of X
function regularize_bcs(B::Array)

    # permute rows of X/columns of B so principle block of B is nonsingular
    P = nonsingular_permute(B)

    B = B*P

    # we apply Q' to upper triangularize B,
    # avoiding any need to permute rows
    Q,R = qr(B)
    Q=Q[:,1:size(B,1)]

    K = size(B,1)

    # we invert the principle block of R
    # so that the BC leads with the identity
    Q = inv(R[:,1:K])*Q'
    R = inv(R[:,1:K])*R

    R,Q,P
end

function regularize_bcs{T<:Number}(B::Array, G::Matrix{T})
    R,Q,P=regularize_bcs(B)
    # we invert the principle block of R
    # so that the BC leads with the identity
    G = Q*G


    R,G,P
end


function regularize_bcs(B::Array, G::Array, L::Array, M::Array)
    if length(B) == 0
        R = B
        P = eye(size(L,2))
    else
        R,G,P=regularize_bcs(B,G)
        L = L*P
        M = M*P
    end

    R,G, L, M, P
end

# Result satisfies
#       inv(Q)*R*P' == B
#       L*P' == L (input)
#       M*P' == M (input)
function regularize_bcs(B::Array, L::Array, M::Array)
    if length(B) == 0
        R = B
        P = eye(size(L,2))
        Q= eye(0)
    else
        R,Q,P=regularize_bcs(B)

        L = L*P
        M = M*P
    end

    R,Q,L,M,P
end


function regularize_bcs(B::Array,Ls::Vector)
   if length(B) == 0
        R = B
        P = eye(size(B,2))
        Q= eye(0)
    else
        R,Q,P=regularize_bcs(B)
        Ls=[L*P for L in Ls]
    end
    R,Q,Ls,P
end


function reduce_dofs( R,G, Mx, My, F )
    if length(R) > 0
        GM = G*My.'
        for k = 1:size(R,1)
            F = F - Mx[:,k]*GM[k:k,:]
            Mx = Mx - Mx[:,k]*R[k:k,:]
        end
    end

    Mx, F
end




function reduce_constrained(Bx,Gx,By,Gy,Lx,Ly,Mx,My,F)
    Rx,Gx,Lx,Mx,Px=regularize_bcs(Bx,Gx,Lx,Mx)
    Ry,Gy,Ly,My,Py=regularize_bcs(By,Gy,Ly,My)

    Lx,F = reduce_dofs(Rx,Gx,Lx,Ly,F)
    Mx,F = reduce_dofs(Rx,Gx,Mx,My,F)
    Ly,F = reduce_dofs(Ry,Gy,Ly,Lx,F.');    F = F.';
    My,F = reduce_dofs(Ry,Gy,My,Mx,F.');    F = F.';

    Kx = size(Bx,1); Ky = size(By,1);
    A=Lx[:,Kx+1:end];B=Ly[:,Ky+1:end];
    C=Mx[:,Kx+1:end];D=My[:,Ky+1:end];

    A,B,C,D,F,Rx,Gx,Px,Ry,Gy,Py
end


constrained_lyap(B,L,M,F)=constrained_lyap(B[1,1],B[1,2],B[2,1],B[2,2],L[1],L[2],M[1],M[2],F)

function constrained_lyap(Bx,Gx,By,Gy,Lx,Ly,Mx,My,F)
    A,B,C,D,F,Rx,Gx,Px,Ry,Gy,Py = reduce_constrained(Bx,Gx,By,Gy,Lx,Ly,Mx,My,F)

    X22=lyap(A,B,C,D,F)

    Kx = size(Bx,1); Ky = size(By,1);
    X12 = Gx[:,Ky+1:end] - Rx[:,Kx+1:end]*X22
    X21 = Gy[:,Kx+1:end].' - X22*Ry[:,Ky+1:end].'
    X11 = Gx[:,1:Ky] - Rx[:,Kx+1:end]*X21
    X11a= Gy[:,1:Kx].' - X12*Ry[:,Ky+1:end].'


    tol = 1e-13;
    bcerr = norm(X11 - X11a);

    if bcerr>tol
       warn("Boundary conditions differ by " * string(bcerr));
    end


    X = [X11 X12; X21 X22];

    X = Px*X*Py.';

end
