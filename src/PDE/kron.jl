
## Represents a list of operators with degrees of freedom reduce
immutable ReducedDiscreteOperators{BT,MT}
    bcP::Array{BT,2}  # permute columns of bcs
    bcQ::Array{BT,2}  # bc normalizer
    bcs::Array{BT,2} # normalized bcs
    
    ops::Vector{SparseMatrixCSC{MT,Int64}}
    opcols::Vector{SparseMatrixCSC{MT,Int64}}
    
    domainspace::FunctionSpace
    rangespace::FunctionSpace
end

numbcs(A::ReducedDiscreteOperators)=size(A.bcs,1)
Base.size(A::ReducedDiscreteOperators,k...)=size(A.ops,k...)

function ReducedDiscreteOperators(Bx,Ls,nx)
    Ls=promotespaces(Ls)
    B,LLs=pdetoarray(Bx,Ls,nx)
    K=length(Bx)
    opcols=SparseMatrixCSC{Float64,Int64}[sparse(L[:,1:K]) for L in LLs]
    R,Q,LLs,P=regularize_bcs(B,LLs)
    ops=SparseMatrixCSC{Float64,Int64}[sparse(cont_reduce_dofs(R,L)[:,K+1:end]) for L in LLs]
    ReducedDiscreteOperators(P,Q,R,ops,opcols,domainspace(Ls[1]),rangespace(Ls[1]))
end

function Base.kron(Ax::ReducedDiscreteOperators,Ay::ReducedDiscreteOperators)
   @assert length(Ax.ops)==length(Ay.ops)

    M=spzeros(Float64,size(Ax.ops[1],1)*size(Ay.ops[1],1),size(Ax.ops[1],2)*size(Ay.ops[1],2))
    for k=1:length(Ax.ops)
        M+=kron(Ax.ops[k],Ay.ops[k])
    end

    M 
end


regularize_bcs(S::ReducedDiscreteOperators,Gy)=length(Gy)==0?Gy:S.bcQ*Gy


function cont_reduce_dofs{T<:Fun}(Ax::ReducedDiscreteOperators,Ay::Vector,G::Vector{T},F)
    G=regularize_bcs(Ax,G)
    cont_reduce_dofs(Ax.opcols,Ay,G,F)
end
function cont_reduce_dofs{T<:Fun}(Ax::ReducedDiscreteOperators,Ay::ReducedDiscreteOperators,G::Vector{T},F)
    G=regularize_bcs(Ax,G)
    cont_reduce_dofs(Ax.opcols,Ay.ops,coefficients(G)[numbcs(Ay)+1:end,:],F)
end


immutable PDEOperatorKron{T}
    opsx::ReducedDiscreteOperators
    opsy::ReducedDiscreteOperators
    indsBx::Vector{Int}
    indsBy::Vector{Int}    
    kron::SparseMatrixCSC{T,Int64}
    domainspace::BivariateFunctionSpace
    rangespace::BivariateFunctionSpace    
end

domainspace(P::PDEOperatorKron,k)=P.domainspace[k]
rangespace(P::PDEOperatorKron,k)=P.rangespace[k]

function PDEOperatorKron(A,nx::Integer,ny::Integer)
    L=A[end]
    indsBx,Bx=findfunctionals(A,1)
    Ax=ReducedDiscreteOperators(Bx,L.ops[:,1],nx)
    indsBy,By=findfunctionals(A,2)
    Ay=ReducedDiscreteOperators(By,L.ops[:,2],ny)
    PDEOperatorKron(Ax,Ay,indsBx,indsBy,kron(Ax,Ay),domainspace(L),rangespace(L))
end

bcinds(A::PDEOperatorKron,k)=k==1?A.indsBx:A.indsBy
numbcs(A::PDEOperatorKron,k)=numbcs(k==1?A.opsx:A.opsy)

Base.kron{T<:PDEOperator}(A::Vector{T},nx::Integer,ny::Integer)=PDEOperatorKron(A,nx,ny)


function pdesolve(A::PDEOperatorKron,f)
    indsBx=bcinds(A,1);indsBy=bcinds(A,2)    
    fx,fy,F=pde_normalize_rhs(A,f)    
end


