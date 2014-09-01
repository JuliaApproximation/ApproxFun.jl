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
    bcQ::Array{T,2}  # bc normalizer
    
    bcs::Array{T,2} # normalized bcs
    
    # C == QRZ',  D == QTZ'
    # where C/D have degrees of freedom removed
    Q::Array{T,2}   
    Z::Array{T,2}
               
    R::Array{T,2}
    T::Array{T,2}
    
    # L[:,1:k] and M[:,1:k]  so we know how the columns are killed
    Lcols::Array{T,2}
    Mcols::Array{T,2}
    
    domainspace::OperatorSpace
    rangespace::OperatorSpace    
end

Base.size(S::OperatorSchur)=size(S.bcP)
Base.size(S::OperatorSchur,k)=size(S.bcP,k)

numbcs(S::OperatorSchur)=size(S.bcQ,1)
domainspace(S::OperatorSchur)=S.domainspace
rangespace(S::OperatorSchur)=S.rangespace
domain(S::OperatorSchur)=domain(domainspace(S))



OperatorSchur{FT<:Functional}(B::Vector{FT},L::UniformScaling,M::Operator,n::Integer)=OperatorSchur(B,ConstantOperator(L),M,n)
OperatorSchur{FT<:Functional}(B::Vector{FT},L::Operator,M::UniformScaling,n::Integer)=OperatorSchur(B,L,ConstantOperator(M),n)
OperatorSchur{FT<:Functional}(B::Vector{FT},L::Operator,M::Operator,n::Integer)=OperatorSchur(pdetoarray(B,L,M,n)...,findmindomainspace([L,M]),findmaxrangespace([L,M]))
OperatorSchur(B,L::SparseMatrixCSC,M::SparseMatrixCSC,ds,rs)=OperatorSchur(B,full(L),full(M),ds,rs)
function OperatorSchur(B::Array,L::Array,M::Array,ds,rs)
    B,Q,L,M,P=regularize_bcs(B,L,M)
    
    K = size(B,1)    
    
    Lcols=L[:,1:K];Mcols=M[:,1:K]    
    
    L=cont_reduce_dofs(B,L)
    M=cont_reduce_dofs(B,M)
    
    C=L[:,K+1:end]
    D=M[:,K+1:end]
    CD=schurfact(C,D)
    Q2=CD[:left];Z2=CD[:right]
    R=CD[:S]; T=CD[:T]
    
    OperatorSchur(P,Q,B,Q2,Z2,R,T, Lcols,Mcols,ds,rs)
end

