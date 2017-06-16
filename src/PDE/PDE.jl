export discretize,timedirichlet

include("KroneckerOperator.jl")

## PDE

lap(d::Space) = Laplacian(d)
lap(d::Domain) = Laplacian(d)
lap(f::Fun) = Laplacian()*f


function Laplacian(d::BivariateSpace,k::Integer)
    Dx2=Derivative(d,[2,0])
    Dy2=Derivative(d,[0,2])
    if k==1
        LaplacianWrapper(Dx2+Dy2,k)
    else
        @assert k > 0
        Δ=Laplacian(d,1)
        LaplacianWrapper(TimesOperator(Laplacian(rangespace(Δ),k-1),Δ),k)
    end
end

Laplacian(d::BivariateDomain,k::Integer) = Laplacian(Space(d),k)
grad(d::ProductDomain) = [Derivative(d,[1,0]),Derivative(d,[0,1])]


function Dirichlet(d::Union{ProductDomain,TensorSpace},k...)
    @assert length(d)==2
    Bx = $op(d[1],k...)
    By = $op(d[2],k...)
    DirichletWrapper(if isempty(Bx)
        I⊗By
    elseif isempty(By)
        Bx⊗I
    else
        [Bx⊗I;I⊗By]
    end, k...)
end


function timedirichlet(d::Union{ProductDomain,TensorSpace})
    @assert length(d.domains)==2
    Bx=Dirichlet(d.domains[1])
    Bt=Dirichlet(d.domains[2])[1]
    [I⊗Bt;Bx⊗I]
end



function *(B::Operator,f::ProductFun)
    if isafunctional(B)
        Fun(space(f,2),map(c->Number(B*c),f.coefficients))
    else
        ProductFun(space(f),map(c->B*c,f.coefficients))
    end
end

*(f::ProductFun,B::Operator)=(B*(f.')).'
