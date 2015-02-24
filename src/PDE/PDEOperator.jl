
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

⊗(A::Vector,B::Operator)=PDEOperator[PDEOperator(Ai,B) for Ai in A]
⊗(A::Operator,B::Vector)=PDEOperator[PDEOperator(A,Bi) for Bi in B]
⊗(A::Vector,B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗(A::UniformScaling,B::Vector)=ConstantOperator(1.0A.λ)⊗B



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

*(B::Functional,f::ProductFun)=Fun(map(c->B*c,f.coefficients),space(f,2))
*(B::BandedOperator,f::ProductFun)=ProductFun(map(c->B*c,f.coefficients),space(f))

*(f::ProductFun,B::Operator)=B*(f.')


function lap(d::Union(ProductDomain,TensorSpace))
    @assert length(d)==2
    Dx=Derivative(d[1])
    Dy=Derivative(d[2])    
    Dx^2⊗I+I⊗Dy^2
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
            [Bx⊗I;I⊗By]
        end
    end
end


function timedirichlet(d::Union(ProductDomain,TensorSpace))
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    Bt=dirichlet(d.domains[2])[1]
    [I⊗Bt;Bx⊗I]
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


isfunctional(B::PDEOperator,k::Integer)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,k])<:Functional
isxfunctional(B::PDEOperator)=isfunctional(B,1)
isyfunctional(B::PDEOperator)=isfunctional(B,2)
ispdeop(B::PDEOperator)=!isxfunctional(B)&&!isyfunctional(B)


function findfunctionals{T<:PDEOperator}(A::Vector{T},k::Integer)
    indsBx=find(f->isfunctional(f,k),A)
    indsBx,Functional[(@assert Ai.ops[1,k==1?2:1]==ConstantOperator{Float64}(1.0); Ai.ops[1,k]) for Ai in A[indsBx]]    
end



    