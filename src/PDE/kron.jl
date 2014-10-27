
## Represents a list of operators with degrees of freedom reduce
type ReducedDiscreteOperators{BT,MT}
    bcP::Array{BT,2}  # permute columns of bcs
    bcQ::Array{BT,2}  # bc normalizer
    bcs::Array{BT,2} # normalized bcs
    
    ops::Vector{SparseMatrixCSC{MT,Int64}}
    opcols::Vector{SparseMatrixCSC{MT,Int64}}
    
    domainspace::FunctionSpace
    rangespace::FunctionSpace
end

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
    cont_reduce_dofs(Ax.opcols,Ay.ops,coefficients(G)[:,length(Ay.bcs)+1:end],F)
end