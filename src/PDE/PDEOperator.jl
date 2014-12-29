
export lap,grad,timedirichlet, ⊗


type PDEOperator
    # 2x2 array 
    # [Lx Ly;
    # [Mx My] 
    # that represents
    # Lx⊗Ly + Mx⊗My
 
    ops::Array{Operator,2}
    domain::BivariateDomain
end
function PDEOperator(LL::Array)
    j=1;  dx=AnyDomain()
    for k=1:size(LL,1)
        d=domain(LL[k,j]) 
        if d != AnyDomain()
            dx=d
            break
        end
    end 
    j=2;  dy=AnyDomain()
    for k=1:size(LL,1)
        d=domain(LL[k,j]) 
        if d != AnyDomain()
            dy=d
            break
        end
    end   
    PDEOperator(LL,dx*dy)
end



PDEOperator(A::Operator,B::Operator)=PDEOperator([A B])
PDEOperator(A::Operator,B::Operator,d::BivariateDomain)=PDEOperator([A B],d)


domainspace(L::PDEOperator,j::Integer)=findmindomainspace(L.ops[:,j])
rangespace(L::PDEOperator,j::Integer)=findmaxrangespace(L.ops[:,j])

for op in (:domainspace,:rangespace)
    @eval begin
        $op(L::PDEOperator)=$op(L,1)⊗$op(L,2)
    end
end


domain(LL::PDEOperator,j::Integer)=domain(LL)[j]
domain(LL::PDEOperator)=LL.domain


Base.transpose(LL::PDEOperator)=PDEOperator(LL.ops[:,end:-1:1])


## Space promotion


function promotedomainspace(P::PDEOperator,S::FunctionSpace,col::Integer)
    ret=copy(P.ops)
    for k=1:size(P.ops,1)
        ret[k,col]=promotedomainspace(ret[k,col],S)
    end
    
    PDEOperator(ret,domain(P))
end
promotedomainspace(P::PDEOperator,S::TensorSpace)=promotedomainspace(promotedomainspace(P,S[1],1),S[2],2)



## Algebra

⊗(A::Operator,B::Operator)=PDEOperator(A,B)
⊗(A::Operator,B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗(A::UniformScaling,B::Operator)=ConstantOperator(1.0A.λ)⊗B

⊗(A::Fun,B::Fun)=Multiplication(A)⊗Multiplication(B)
⊗(A,B::Fun)=A⊗Multiplication(B)
⊗(A::Fun,B)=Multiplication(A)⊗B

⊗{T<:Operator}(A::Vector{T},B::Operator)=PDEOperator[PDEOperator(Ai,B) for Ai in A]
⊗{T<:Operator}(A::Operator,B::Vector{T})=PDEOperator[PDEOperator(A,Bi) for Bi in B]
⊗{T<:Operator}(A::Vector{T},B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗{T<:Operator}(A::UniformScaling,B::Vector{T})=ConstantOperator(1.0A.λ)⊗B



function +(A::PDEOperator,B::PDEOperator)

    ret=copy(A.ops)
    for k=1:size(B.ops,1)
        ##TODO: might not be ordered
        if ret[k,1]==B.ops[k,1]
            ret[k,2]+=B.ops[k,2]
        elseif ret[k,2]==B.ops[k,2]
            ret[k,1]+=B.ops[k,1]            
        else
            ret=[ret;B.ops[k,:]]
        end
    end
    
    # the following is so operators know
    # they are on Disk
    if !isa(domain(A),ProductDomain)
        PDEOperator(ret,domain(A))
    elseif !isa(domain(B),ProductDomain)
        PDEOperator(ret,domain(B))
    else
        # both are product domains, reconstruct the domain
        PDEOperator(ret) 
    end
end


function lap(d::Union(ProductDomain,TensorSpace))
    @assert length(d)==2
    Dx=Derivative(d[1])
    Dy=Derivative(d[2])    
    Dx^2⊗I+I⊗Dy^2
end



function -(A::PDEOperator)
    ops=copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=-ops[k,1]
    end
    PDEOperator(ops,domain(A))
end

+(A::UniformScaling,B::PDEOperator)=B+PDEOperator(ConstantOperator(1.0A.λ),ConstantOperator(1.0),domain(B))
+(B::PDEOperator,A::UniformScaling)=B+PDEOperator(ConstantOperator(1.0A.λ),ConstantOperator(1.0),domain(B))
-(A::UniformScaling,B::PDEOperator)=-B+PDEOperator(ConstantOperator(1.0A.λ),ConstantOperator(1.0),domain(B))
-(B::PDEOperator,A::UniformScaling)=B+PDEOperator(ConstantOperator(-1.0A.λ),ConstantOperator(1.0),domain(B))

-(A::PDEOperator,B::PDEOperator)=A+(-B)

function *(c::Number,A::PDEOperator)
    ops=copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=c*ops[k,1]
    end
    PDEOperator(ops,domain(A))
end
*(A::PDEOperator,c::Number)=c*A
function *(A::PDEOperator,B::PDEOperator)
    # TODO: higher rank operators
    if size(A.ops,1)==size(B.ops,1)==1
        @assert size(A.ops,2)==size(B.ops,2)==2
        ops=[A.ops[1,1]*B.ops[1,1] A.ops[1,2]*B.ops[1,2]]
    else
        @assert A==B
        @assert size(A.ops,1)==size(B.ops,1)==2
        @assert size(A.ops,2)==size(B.ops,2)==2 
        
        ops=[A.ops[1,1]^2 A.ops[1,2]^2;
                     A.ops[1,1]*A.ops[2,1] A.ops[1,2]*A.ops[2,2];
                     A.ops[2,1]*A.ops[1,1] A.ops[2,2]*A.ops[1,2];                     
                     A.ops[2,1]^2 A.ops[2,2]^2
                    ]
    end
    
    if !isa(domain(A),ProductDomain)
        PDEOperator(ops,domain(A))
    elseif !isa(domain(B),ProductDomain)
        PDEOperator(ops,domain(B))
    else
        # both are product domains, reconstruct the domain
        PDEOperator(ops) 
    end
end



##TODO how to determine whether x or y?
function *(a::Fun,A::PDEOperator)
    ops = copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=a*ops[k,1]
    end
    PDEOperator(ops,domain(A))
end

Derivative(d::Union(ProductDomain,TensorSpace),k::Integer)=k==1?Derivative(d[1])⊗I:I⊗Derivative(d[2])  
Base.diff(d::Union(ProductDomain,TensorSpace),k::Integer)=Derivative(d,k)
grad(d::ProductDomain)=[Derivative(d,k) for k=1:length(d.domains)]



for op in (:dirichlet,:neumann,:diffbcs)
    @eval begin
        function $op(d::Union(ProductDomain,TensorSpace),k...)
            @assert length(d)==2
            Bx=$op(d[1],k...)
            By=$op(d[2],k...)
            [Bx⊗I,I⊗By]
        end
    end
end


function timedirichlet(d::Union(ProductDomain,TensorSpace))
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    Bt=dirichlet(d.domains[2])[1]
    [I⊗Bt,Bx⊗I]
end



function *{S,T}(L::PDEOperator,f::LowRankFun{S,T})
    @assert size(L.ops,2)==2
    @assert size(L.ops,1)==2    
    n=length(f.A)
    A=Array(Fun,2n)
    B=Array(Fun,2n)
    
    for k=1:n
        A[k]=L.ops[1,1]*f.A[k]
        B[k]=L.ops[1,2]*f.B[k]        
        A[k+n]=L.ops[2,1]*f.A[k]
        B[k+n]=L.ops[2,2]*f.B[k]        
    end
    
    LowRankFun(A,B)
end

function tensormult(A::BandedOperator,B::BandedOperator,F::ProductFun)
    ##TODO: general bands?
    B=promotedomainspace(B,space(F,2))
    @assert bandinds(B)==(0,0)
    @assert rangespace(B)==domainspace(B)
    ret=copy(F.coefficients)
    for k=1:length(ret)
        ret[k]=B[k,k]*A*ret[k]
    end
    ProductFun(ret,space(F))
end


function *(A::PDEOperator,F::ProductFun)
    ret=tensormult(A.ops[1,1],A.ops[1,2],F)
    for k=2:size(A.ops,1)
        ret+=tensormult(A.ops[k,1],A.ops[k,2],F)
    end
    ret
end




## Schur PDEOperator
# represent an operator that's been discretized in the Y direction

abstract AbstractPDEOperatorSchur

immutable PDEOperatorSchur{OS<:AbstractOperatorSchur,LT<:Number,MT<:Number,ST<:Number,FT<:Functional} <: AbstractPDEOperatorSchur
    Bx::Vector{FT}
    Lx::Operator{LT}
    Mx::Operator{MT}
    S::OS
    
    indsBx::Vector{Int}
    indsBy::Vector{Int}
    
    Rdiags::Vector{SavedBandedOperator{ST}}
end

#TODO: Functional type?
Base.eltype{OS,LT,MT,ST,FT}(PO::PDEOperatorSchur{OS,LT,MT,ST,FT})=promote_type(eltype(PO.S),LT,MT,ST)


function PDEOperatorSchur{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy)
    Xops=promotespaces([Lx,Mx])
    Lx=SavedBandedOperator(Xops[1]);Mx=SavedBandedOperator(Xops[2])

    ny=size(S,1)
    nbcs=numbcs(S)
    Rdiags=Array(SavedBandedOperator{promote_type(LT,MT,BT,ST)},ny)
    
    resizedata!(Lx,ny);resizedata!(Mx,ny)
    
    
    for k=1:ny-nbcs
        ##TODO: Do block case
        Rdiags[k]=SavedBandedOperator(getdiagonal(S,k,1)*Lx + getdiagonal(S,k,2)*Mx)
        resizedata!(Rdiags[k],ny)
    end
    
    if isempty(Bx)
        funcs=Array(SavedFunctional{Float64},0)
    else
        funcs=Array(SavedFunctional{mapreduce(eltype,promote_type,Bx)},length(Bx))
        for k=1:length(Bx)
            funcs[k]=SavedFunctional(Bx[k])
            resizedata!(funcs[k],ny)
        end
    end

    
    PDEOperatorSchur(funcs,Lx,Mx,S,indsBx,indsBy,Rdiags)
end




PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::Operator,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,Lx,Mx,S,[1:length(Bx)],length(Bx)+[1:numbcs(S)])
PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::UniformScaling,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,Lx,ConstantOperator(Mx.λ),S)
PDEOperatorSchur(Bx::Vector,Lx::UniformScaling,Mx::Operator,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,ConstantOperator(Lx.λ),Mx,S)


function PDEOperatorSchur(Bx,By,A::PDEOperator,ny::Integer,indsBx,indsBy)
   @assert size(A.ops)==(2,2)
   PDEOperatorSchur(Bx,A.ops[1,1],A.ops[2,1],schurfact(By,A.ops[:,2],ny),indsBx,indsBy)
end

isfunctional(B::PDEOperator,k::Integer)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,k])<:Functional
isxfunctional(B::PDEOperator)=isfunctional(B,1)
isyfunctional(B::PDEOperator)=isfunctional(B,2)
ispdeop(B::PDEOperator)=!isxfunctional(B)&&!isyfunctional(B)


bcinds(S::PDEOperatorSchur,k)=k==1?S.indsBx:S.indsBy
numbcs(S::AbstractPDEOperatorSchur,k)=length(bcinds(S,k))
Base.length(S::PDEOperatorSchur)=size(S.S,1)

function findfunctionals{T<:PDEOperator}(A::Vector{T},k::Integer)
    indsBx=find(f->isfunctional(f,k),A)
    indsBx,Functional[(@assert Ai.ops[1,k==1?2:1]==ConstantOperator{Float64}(1.0); Ai.ops[1,k]) for Ai in A[indsBx]]    
end



function PDEOperatorSchur{T<:PDEOperator}(A::Vector{T},ny::Integer)
    indsBx,Bx=findfunctionals(A,1)
    indsBy,By=findfunctionals(A,2)    

    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)    
    LL=A[end]
    
    
    PDEOperatorSchur(Bx,By,LL,ny,indsBx,indsBy)
end



for op in (:domainspace,:rangespace)
    @eval begin
        $op(P::PDEOperatorSchur,k::Integer)=k==1?$op(P.Lx):$op(P.S)
        $op(L::PDEOperatorSchur)=$op(L,1)⊗$op(L,2)
    end
end


domain(P::PDEOperatorSchur,k::Integer)=k==1?domain(P.Lx):domain(P.S)
domain(P::PDEOperatorSchur)=domain(domainspace(P))



#############
# PDEStrideOperatorSchur
#############


immutable PDEStrideOperatorSchur  <: AbstractPDEOperatorSchur
    odd    
    even
end


function PDEStrideOperatorSchur(Bx,Lx::Operator,Mx::Operator,S::StrideOperatorSchur,indsBx,indsBy)
    @assert indsBy==[3,4]
    Podd=PDEOperatorSchur(Bx,Lx,Mx,S.odd,indsBx,[indsBy[1]])
    Peven=PDEOperatorSchur(Bx,Lx,Mx,S.even,indsBx,[indsBy[1]])
    
    PDEStrideOperatorSchur(Podd,Peven)
end

## Product

# Represents an operator on e.g. a Disk
immutable PDEProductOperatorSchur{ST<:Number,FT<:Functional} <: AbstractPDEOperatorSchur
    Bx::Vector{Vector{FT}}
    Rdiags::Vector{SavedBandedOperator{ST}}
    domainspace::AbstractProductSpace
    indsBx::Vector{Int}
end

Base.length(S::PDEProductOperatorSchur)=length(S.Rdiags)

function PDEProductOperatorSchur{T<:PDEOperator}(A::Vector{T},sp::AbstractProductSpace,nt::Integer)
    indsBx=find(isxfunctional,A)
    Bx=Functional[(@assert Ai.ops[1,2]==ConstantOperator{Float64}(1.0); Ai.ops[1,1]) for Ai in A[indsBx]]
    
    
    indsBy=find(isyfunctional,A)
    @assert length(indsBy)==0
    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)
    
    L=A[end]
    
    
    
    L=promotedomainspace(L,sp[2],2)
    S=schurfact([],L.ops[:,2],nt)
    @assert(isa(S,DiagonalOperatorSchur))

    #TODO: Space Type
    BxV=Array(Vector{SavedFunctional{Float64}},nt)
    Rdiags=Array(SavedBandedOperator{Complex{Float64}},nt)    
    for k=1:nt
        op=getdiagonal(S,k,1)*L.ops[1,1]
        for j=2:size(L.ops,1)
            Skj=getdiagonal(S,k,j)
            if Skj!=0   # don't add unnecessary ops
                op+=Skj*L.ops[j,1]         
            end
        end
    
        csp=columnspace(sp,k)
        Rdiags[k]=SavedBandedOperator(promotedomainspace(op,csp))
        BxV[k]=SavedFunctional{Float64}[SavedFunctional(promotedomainspace(Bxx,csp)) for Bxx in Bx]
        resizedata!(Rdiags[k],2nt+100)
        for Bxx in BxV[k]
            resizedata!(Bxx,2nt+100)        
        end
    end  
    PDEProductOperatorSchur(BxV,Rdiags,sp,indsBx)
end

PDEProductOperatorSchur{T<:PDEOperator}(A::Vector{T},sp::BivariateDomain,nt::Integer)=PDEProductOperatorSchur(A,Space(sp),nt)


##TODO: Larger Bx
bcinds(S::PDEProductOperatorSchur,k)=k==1?S.indsBx:[]
domainspace(S::PDEProductOperatorSchur)=S.domainspace
domainspace(S::PDEProductOperatorSchur,k)=S.domainspace[k]

# for op in (:domainspace,:rangespace)
#     @eval begin
#         $op(P::PDEProductOperatorSchur,k::Integer)=k==1?$op(P.Lx):$op(P.S)
#        $op(L::PDEProductOperatorSchur)=$op(L,1)⊗$op(L,2)
#     end
# end






type ProductRangeSpace <: BivariateFunctionSpace
    S::PDEProductOperatorSchur
end

rangespace(S::PDEProductOperatorSchur)=ProductRangeSpace(S)

function space(S::ProductRangeSpace,k)
    @assert k==2
    S.S.domainspace[2]
end 

function coefficients{S,V,SS,T}(f::ProductFun{S,V,SS,T},sp::ProductRangeSpace)
    @assert space(f,2)==space(sp,2)
    
    n=min(size(f,2),length(sp.S))
    F=[coefficients(f.coefficients[k],rangespace(sp.S.Rdiags[k])) for k=1:n]
    m=mapreduce(length,max,F)
    ret=zeros(T,m,n)
    for k=1:n
        ret[1:length(F[k]),k]=F[k]
    end
    ret    
end



## Constructuor


Base.schurfact{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy)=PDEOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy)
Base.schurfact{LT<:Number,MT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::StrideOperatorSchur{ST},indsBx,indsBy)=PDEStrideOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy)

function Base.schurfact(Bx,By,A::PDEOperator,ny::Integer,indsBx,indsBy)
   @assert size(A.ops)==(2,2)
   schurfact(Bx,A.ops[1,1],A.ops[2,1],schurfact(By,A.ops[:,2],ny),indsBx,indsBy)
end



function Base.schurfact{T<:PDEOperator}(A::Vector{T},S::ProductDomain,ny::Integer)
    indsBx,Bx=findfunctionals(A,1)
    indsBy,By=findfunctionals(A,2)    

    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)    
    LL=A[end]
    
    
    schurfact(Bx,By,LL,ny,indsBx,indsBy)
end





Base.schurfact{T<:PDEOperator}(A::Vector{T},S::BivariateDomain,n::Integer)=PDEProductOperatorSchur(A,S,n)


Base.schurfact{T<:PDEOperator}(A::Vector{T},n::Integer)=schurfact(A,domain(A[end]),n)
Base.schurfact(A::PDEOperator,n::Integer)=schurfact([A],n)

function *(A::PDEProductOperatorSchur,F::ProductFun)
    ret=copy(F.coefficients)
    for k=1:length(ret)
        ret[k]=A.Rdiags[k]*ret[k]
    end
    ProductFun(ret,space(F))
end




## discretize

#TODO: don't hard code Disk
discretize{T<:PDEOperator}(A::Vector{T},S...)=size(A[end].ops,1)==2||domain(A[end])==Disk()?schurfact(A,S...):kron(A,S...)
discretize(A::PDEOperator,n::Integer)=discretize([A],n)


    