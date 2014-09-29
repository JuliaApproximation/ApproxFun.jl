##TODO: Unify with coefficients()
toarray{T<:Functional}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
toarray{T<:Number}(B::Array{Fun{T}},n)=T[    j<=length(B[k])?B[k].coefficients[j]:0 for  k=1:length(B),j=1:n]

iscomplexfunornumber(A)=false
iscomplexfunornumber(A::Complex{Float64})=true
iscomplexfunornumber{S}(A::Fun{S,Complex{Float64}})=true

function toarray(B::Array,n)
    T=Float64
    for Bk in B
        if iscomplexfunornumber(Bk)
            T=Complex{Float64}
        end
    end    

    ret = zeros(T,length(B),n)
    
    for k=1:length(B), j=1:n
        if  isa(B[k],Fun)
            ret[k,j] = j<=length(B[k])?B[k].coefficients[j]:0 
        elseif isa(B[k],Number) && j == 1
            ret[k,j] = B[k]
        end
    end

    ret
end

function toarray{T<:Operator}(A::Vector{T},n::Integer,m::Integer)
    ret = zeros(n,m)
    
    nbc = isa(A[end],Functional)?length(A):length(A)-1
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
        Q= eye(0)
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


abstract AbstractOperatorSchur{BT,MT}

domainspace(S::AbstractOperatorSchur)=S.domainspace
rangespace(S::AbstractOperatorSchur)=S.rangespace
domain(S::AbstractOperatorSchur)=domain(domainspace(S))


type DiagonalOperatorSchur{MT<:Number} <:AbstractOperatorSchur{MT,MT}
    R::Vector{MT}
    T::Vector{MT}
    
    domainspace::FunctionSpace
    rangespace::FunctionSpace       
end

function DiagonalOperatorSchur{T1<:Number,T2<:Number}(R::Vector{T1},T::Vector{T2},d,r)
    PT=promote_type(T1,T2)
    DiagonalOperatorSchur(convert(Vector{PT},R),convert(Vector{PT},T),d,r)
end

function DiagonalOperatorSchur(L::BandedOperator,M::BandedOperator,n::Integer)
    Yop=promotespaces([L,M])    
    DiagonalOperatorSchur(diag(Yop[1][1:n,1:n]),diag(Yop[2][1:n,1:n]),domainspace(Yop[1]),rangespace(Yop[2]))
end

Base.size(S::DiagonalOperatorSchur,k)=length(S.R)
Base.size(S::DiagonalOperatorSchur)=size(S,1),size(S,2)

getdiagonal(S::DiagonalOperatorSchur,k,j)=j==1?S.R[k]:S.T[k]

numbcs(::DiagonalOperatorSchur)=0


type OperatorSchur{BT<:Number,MT<:Number} <:AbstractOperatorSchur{BT,MT}
    bcP::Array{BT,2}  # permute columns of bcs
    bcQ::Array{BT,2}  # bc normalizer
    
    bcs::Array{BT,2} # normalized bcs
    
    # C == QRZ',  D == QTZ'
    # where C/D have degrees of freedom removed
    Q::Array{MT,2}   
    Z::Array{MT,2}
               
    R::Array{MT,2}
    T::Array{MT,2}
    
    # L[:,1:k] and M[:,1:k]  so we know how the columns are killed
    Lcols::Array{MT,2}
    Mcols::Array{MT,2}
    
    domainspace::FunctionSpace
    rangespace::FunctionSpace    
end

#make sure cols are same type as ops
OperatorSchur{M<:Number}(bcP,bcQ,bcs,Q::Array{M,2},Z::Array{M,2},R::Array{M,2},T::Array{M,2},Lcols::Array,Mcols::Array,ds,rs)=
    OperatorSchur(bcP,bcQ,bcs,Q,Z,R,T,convert(Array{M},Lcols),convert(Array{M},Mcols),ds,rs)

Base.size(S::OperatorSchur)=size(S.bcP)
Base.size(S::OperatorSchur,k)=size(S.bcP,k)

numbcs(S::OperatorSchur)=size(S.bcQ,1)

getdiagonal(S::OperatorSchur,k,j)=j==1?S.R[k,k]:S.T[k,k]



OperatorSchur{FT<:Functional}(B::Vector{FT},L::UniformScaling,M::Operator,n::Integer)=OperatorSchur(B,ConstantOperator(L),M,n)
OperatorSchur{FT<:Functional}(B::Vector{FT},L::Operator,M::UniformScaling,n::Integer)=OperatorSchur(B,L,ConstantOperator(M),n)
function OperatorSchur{FT<:Functional}(B::Vector{FT},L::Operator,M::Operator,n::Integer)
    if isempty(B) && bandinds(L)==bandinds(M)==(0,0)
        DiagonalOperatorSchur(L,M,n)
    else
        OperatorSchur(pdetoarray(B,L,M,n)...,findmindomainspace([L,M]),findmaxrangespace([L,M]))
    end
end
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

