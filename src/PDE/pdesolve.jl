export pdesolve

include("lyap.jl")


toarray{T<:RowOperator}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
toarray{T<:IFun}(B::Array{T},n)=Float64[    j<=length(B[k])?B[k].coefficients[j]:0 for  k=1:length(B),j=1:n];
function toarray(B::Any,n)
    ret = zeros(length(B),n)
    
    for k=1:length(B), j=1:n
        if  typeof(B[k]) <: IFun
            ret[k,j] = j<=length(B[k])?B[k].coefficients[j]:0 
        elseif typeof(B[k]) <: Number && j == 1
            ret[k,j] = B[k]
        end
    end

    ret
end


function nonsingular_permute(B)
    K = size(B,1)
    
    k = 1
    
    while rank(B[:,k:K+k-1]) < K
        k= k+1
    end
    
    P = eye(size(B,2))
    
    P = P[:,[k:K+k-1,1:k-1,K+k:end]]
end


function regularize_bcs(B, G, L, M)
    if length(B) == 0
        R = B
        P = eye(size(L,2))
    else
        P = nonsingular_permute(B)
        
        B = B*P
        
        L = L*P
        M = M*P
        
        Q,R = qr(B)
        Q=Q[:,1:size(B,1)]
        
        K = size(B,1)
        
        G = inv(R[:,1:K])*Q'*G
        R = inv(R[:,1:K])*R
    end
    
    R,G, L, M, P
end


function reduce_dofs( R,G, Mx, My, F )
    if length(R) > 0
        GM = G*My.'
        for k = 1:size(R,1)
                F = F - Mx[:,k]*GM[k,:]
                Mx = Mx - Mx[:,k]*R[k,:]
        end
    end
        
    Mx, F
end







constrained_lyap(B,L,M,F)=constrained_lyap(B[1,1],B[1,2],B[2,1],B[2,2],L[1],L[2],M[1],M[2],F)

function constrained_lyap(Bx,Gx,By,Gy,Lx,Ly,Mx,My,F)
     Rx,Gx,Lx,Mx,Px=regularize_bcs(Bx,Gx,Lx,Mx)
    Ry,Gy,Ly,My,Py=regularize_bcs(By,Gy,Ly,My)

    Lx,F = reduce_dofs(Rx,Gx,Lx,Ly,F)
    Mx,F = reduce_dofs(Rx,Gx,Mx,My,F)
    Ly,F = reduce_dofs(Ry,Gy,Ly,Lx,F.');    F = F.';
    My,F = reduce_dofs(Ry,Gy,My,Mx,F.');    F = F.';       

    Kx = size(Bx,1); Ky = size(By,1);
    A=Lx[:,Kx+1:end];B=Ly[:,Ky+1:end];
    C=Mx[:,Kx+1:end];D=My[:,Ky+1:end];

    X22=lyap(A,B,C,D,F)


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

pdesolve(Bxin,Byin,Lin,Min,F,n::Integer)=pdesolve(Bxin,Byin,Lin,Min,F,n,n)

function pdesolve(Bxin,Byin,Lin,Min,Fin::Number,nx::Integer,ny::Integer)
    F=zeros(nx-length(Bxin[1]),ny-length(Byin[1]))
    F[1,1]=Fin
    
    pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
end

function pdesolve(Bxin,Byin,Lin,Min,Fin::Fun2D,nx::Integer,ny::Integer)
    Xop=promotespaces([Lin[1],Min[1]])
    Yop=promotespaces([Lin[2],Min[2]])
    
    Xsp=rangespace(Xop[1])
    Ysp=rangespace(Yop[1])    
    
    nbcx=length(Bxin[1]);nbcy=length(Byin[1])    
    
    F=pad(coefficients(Fin,Xsp,Ysp),nx-nbcx,ny-nbcy)
    pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
end

function pdesolve(Bxin,Byin,Lin,Min,F::Array,nx::Integer,ny::Integer)

    Xop=promotespaces([Lin[1],Min[1]])
    Yop=promotespaces([Lin[2],Min[2]])

    
    Bx=toarray(Bxin[1],nx);By=toarray(Byin[1],ny)
    Gx=toarray(Bxin[2],ny);Gy=toarray(Byin[2],nx)    
    
    nbcx=length(Bxin[1]);nbcy=length(Byin[1]);
    
    Lx=Xop[1][1:nx-nbcx,1:nx];Ly=Yop[1][1:ny-nbcy,1:ny]
    Mx=Xop[2][1:nx-nbcx,1:nx];My=Yop[2][1:ny-nbcy,1:ny]    
    
    
    Fun2D(constrained_lyap({Bx Gx; By Gy},{Lx,Ly},{Mx,My},F),
            domain(Lin[1]),domain(Lin[2]))
end

