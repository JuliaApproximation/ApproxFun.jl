export discretize,timedirichlet


# Bivariate functions have BandedMatrix
op_eltype{T,D}(sp::Space{T,D,2})=BandedMatrix{promote_type(eltype(sp),eltype(domain(sp)))}
op_eltype_realdomain{T,D}(sp::Space{T,D,2})=BandedMatrix{promote_type(eltype(sp),real(eltype(domain(sp))))}

include("OperatorSchur.jl")
include("KroneckerOperator.jl")
include("dekron.jl")
include("diagop.jl")

include("factorizations.jl")
include("pdesolve.jl")
include("kron.jl")

## PDE

lap(d)=Laplacian(d)


function Laplacian(d::Union{ProductDomain,TensorSpace},k::Integer)
    @assert length(d)==2
    Dx2=Derivative(d,[2,0])
    Dy2=Derivative(d,[0,2])
    if k==1
        LaplacianWrapper(Dx2+Dy2,k)
    else
        LaplacianWrapper((Dx2+Dy2)^k,k)
    end
end


grad(d::ProductDomain)=[Derivative(d,[1,0]),Derivative(d,[0,1])]



for op in (:dirichlet,:neumann,:diffbcs)
    @eval begin
        function $op(d::Union{ProductDomain,TensorSpace},k...)
            @assert length(d)==2
            Bx=$op(d[1],k...)
            By=$op(d[2],k...)
            if isempty(Bx)
                I⊗By
            elseif isempty(By)
                Bx⊗I
            else
                [Bx⊗I;I⊗By]
            end
        end
    end
end


function timedirichlet(d::Union{ProductDomain,TensorSpace})
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    Bt=dirichlet(d.domains[2])[1]
    [I⊗Bt;Bx⊗I]
end



function *(B::Operator,f::ProductFun)
    if isafunctional(B)
        Fun(map(c->B*c,f.coefficients),space(f,2))
    else
        ProductFun(map(c->B*c,f.coefficients),space(f))
    end
end

*(f::ProductFun,B::Operator)=(B*(f.')).'
