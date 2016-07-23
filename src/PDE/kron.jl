export kronfact


## Represents a list of operators with degrees of freedom reduce
immutable ReducedDiscreteOperators{BT,MT}
    bcP::Array{BT,2}  # permute columns of bcs
    bcQ::Array{BT,2}  # bc normalizer
    bcs::Array{BT,2} # normalized bcs

    ops::Vector{SparseMatrixCSC{MT,Int64}}
    opcols::Vector{SparseMatrixCSC{MT,Int64}}

    domainspace::Space
    rangespace::Space
end

numbcs(A::ReducedDiscreteOperators)=size(A.bcs,1)
Base.size(A::ReducedDiscreteOperators,k...)=size(A.ops[1],k...)

function ReducedDiscreteOperators(Bx,Ls,nx)
    Ls=promotespaces(Ls)
    B,LLs=pdetoarray(Bx,Ls,nx)
    K=length(Bx)
    opcols=SparseMatrixCSC{Float64,Int64}[sparse(L[:,1:K]) for L in LLs]
    R,Q,LLs,P=regularize_bcs(B,LLs)
    ops=SparseMatrixCSC{Float64,Int64}[sparse(cont_reduce_dofs(R,full(L))[:,K+1:end]) for L in LLs]
    ReducedDiscreteOperators(P,Q,R,ops,opcols,domainspace(Ls[1]),rangespace(Ls[1]))
end

function Base.kron(Ax::ReducedDiscreteOperators,Ay::ReducedDiscreteOperators)
   @assert length(Ax.ops)==length(Ay.ops)

    M=spzeros(Float64,size(Ax,1)*size(Ay,1),size(Ax,2)*size(Ay,2))
    for k=1:length(Ax.ops)
        M+=kron(Ax.ops[k],Ay.ops[k])
    end

    M
end


regularize_bcs(S::ReducedDiscreteOperators,Gy)=length(Gy)==0?Gy:S.bcQ*Gy


function cont_reduce_dofs!{T<:Fun}(Ax::ReducedDiscreteOperators,Ay::Vector,G::Vector{T},F)
    G=regularize_bcs(Ax,G)
    cont_reduce_dofs!(Ax.opcols,Ay,G,F)
end
function cont_reduce_dofs!{T<:Fun}(Ax::ReducedDiscreteOperators,Ay::ReducedDiscreteOperators,G::Vector{T},F)
    G=regularize_bcs(Ax,G)
    cont_reduce_dofs!(Ax.opcols,Ay.ops,coefficients(G)[numbcs(Ay)+1:end,:],F)
end


immutable PDEOperatorKron{T}
    op::Operator{T} #we need to remember the op to reduce DOFs
    opsx::ReducedDiscreteOperators
    opsy::ReducedDiscreteOperators
    indsBx::Vector{Int}
    indsBy::Vector{Int}
    kron::SparseMatrixCSC{T,Int64}
end

Base.eltype{T}(::PDEOperatorKron{T})=T

for op in (:domainspace,:rangespace)
    @eval $op(P::PDEOperatorKron,k...)=$op(P.op,k...)
end

domain(P::PDEOperatorKron)=domain(domainspace(P))

function PDEOperatorKron(A,nx::Integer,ny::Integer)
    L=A[end]
    @assert iskronsumop(L)
    ops=sumops(L)

    indsBx,Bx=findfunctionals(A,1)
    Ax=ReducedDiscreteOperators(Bx,map(op->dekron(op,1),ops),nx)
    indsBy,By=findfunctionals(A,2)
    Ay=ReducedDiscreteOperators(By,map(op->dekron(op,2),ops),ny)
    PDEOperatorKron(L,Ax,Ay,indsBx,indsBy,kron(Ax,Ay))
end

#TODO: could be confused with size(A.kron)...
Base.size(A::PDEOperatorKron,k)=k==1?size(A.opsx,1):size(A.opsy,1)
Base.size(A::PDEOperatorKron)=size(A,1),size(A,2)
bcinds(A::PDEOperatorKron,k)=k==1?A.indsBx:A.indsBy
numbcs(A::PDEOperatorKron,k)=numbcs(k==1?A.opsx:A.opsy)

kronfact(A::Vector,nx::Integer,ny::Integer)=PDEOperatorKron(A,nx,ny)
kronfact(A::Vector,n::Integer)=kronfact(A,n,n)
kronfact(A::Vector,S::BivariateDomain,n::Integer)=kronfact(A,n)


function pdesolve(K::PDEOperatorKron,G)
    fx,fy,F=pde_standardize_rhs(K,G)

    F=cont_reduce_dofs!(K.opsx,map(op->dekron(op,2),sumops(K.op)),fx,F.').'
    F=cont_reduce_dofs!(K.opsy,K.opsx,fy,F)



    X22=reshape(K.kron\vec(pad(coefficients(F,rangespace(K)),size(K)...)),size(K)...)
    Kx,Ky=numbcs(K,1),numbcs(K,2) # doesn't include boundary rows...could be confusing
    nx,ny=size(K,1)+Kx,size(K,2)+Ky

    Gx=pad(coefficients(regularize_bcs(K.opsx,fx)).',:,ny)
    Gy=pad(coefficients(regularize_bcs(K.opsy,fy)).',:,nx)
    Rx,Ry=K.opsx.bcs,K.opsy.bcs

    Px,Py=K.opsx.bcP,K.opsy.bcP

    # following is copyed from constrained_lyap
    X12=Gx[:,Ky+1:end]-Rx[:,Kx+1:end]*X22
    X21=Gy[:,Kx+1:end].'-X22*Ry[:,Ky+1:end].'
    X11 = Gx[:,1:Ky] - Rx[:,Kx+1:end]*X21
    X11a= Gy[:,1:Kx].' - X12*Ry[:,Ky+1:end].'

    tol = 1e-13
    if !isempty(X11)
        bcerr = norm(X11 - X11a)

        if bcerr>tol
           warn("Boundary conditions differ by " * string(bcerr))
        end
    end

    X = [X11 X12; X21 X22]

    ProductFun(Px*X*Py.',domainspace(K))
end


\(A::PDEOperatorKron,G)=pdesolve(A,G)
