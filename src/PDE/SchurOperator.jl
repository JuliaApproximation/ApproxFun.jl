##TODO: Unify with coefficients()
toarray{T<:Functional}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
toarray{T<:IFun}(B::Array{T},n)=Float64[    j<=length(B[k])?B[k].coefficients[j]:0 for  k=1:length(B),j=1:n];
function toarray(B::Array,n)
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

function toarray{T<:Operator}(A::Vector{T},n::Integer,m::Integer)
    ret = zeros(n,m)
    
    nbc = typeof(A[end])<:Functional?length(A):length(A)-1
    for k=1:nbc
        ret[k,:]=A[k][1:m]
    end
    
    if nbc < length(A)
        ret[nbc+1:end,:]=A[end][1:n-nbc,1:m]
    end
    
    ret
end

function pdetoarray(Byin,Lin,Min,ny::Integer)
    Yop=promotespaces([Lin,Min])

    
    By=toarray(Byin,ny)   
    nbcy=length(Byin)
    
    Ly=Yop[1][1:ny-nbcy,1:ny]
    My=Yop[2][1:ny-nbcy,1:ny]    
    
    By,Ly,My    
end

# Result satisfies
#       inv(Q)*R*P' == B
#       L*P' == L (input)
#       M*P' == M (input)
function regularize_bcs(B::Array, L::Array, M::Array)
    if length(B) == 0
        R = B
        P = eye(size(L,2))
    else
        # permute rows of X/columns of B so principle block of B is nonsingular
        P = nonsingular_permute(B)
        
        B = B*P
        
        L = L*P
        M = M*P
        
        # we apply Q' to upper triangularize B,
        # avoiding any need to permute rows
        Q,R = qr(B)
        Q=Q[:,1:size(B,1)]
        
        K = size(B,1)
        
        # we invert the principle block of R
        # so that the BC leads with the identity
        Q = inv(R[:,1:K])*Q'
        R = inv(R[:,1:K])*R
    end
    
    R,Q,L,M,P
end

function cont_reduce_dofs( R,A::Array )
    if length(R) > 0        
        for k = 1:size(R,1)
            A = A - A[:,k]*R[k,:]
        end
    end
        
    A
end


type OperatorSchur{T}
    bcP::Array{T,2}  # permute columns of bcs
    bcQ::Array{T,2}  # normalize bcs
    
    # A == QRZ',  B == QTZ'
    # where A/B have degrees of freedom removed
    Q::Array{T,2}   
    Z::Array{T,2}
               
    R::Array{T,2}
    T::Array{T,2}
end


function OperatorSchur{FT<:Functional}(B::Vector{FT},L::Operator,M::Operator,n::Integer)
    B,L,M=pdetoarray(B,L,M,n)
    R,Q,L,M,P=regularize_bcs(B,full(L),full(M))
    
    
    L=cont_reduce_dofs(R,L)
    M=cont_reduce_dofs(R,M)
    
    K = size(B,1)
    B=L[:,K+1:end]
    D=M[:,K+1:end]
    BD=schurfact(B,D)
    Q2=BD[:left];Z2=BD[:right]
    R=BD[:S]; T=BD[:T]
    
    OperatorSchur(P,Q,Q2,Z2,R,T)
end

