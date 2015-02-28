#####
# This includes short-hands that are convenient but perhaps confusing
#####


export ∫,⨜,⨍

# diff(::Fun{ArraySpace}) is left as array diff

Base.diff{S,T}(f::Fun{S,T},n...)=differentiate(f,n...)
Base.diff(u::MultivariateFun,j...)=differentiate(u,j...)

Base.diff(d::FunctionSpace,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain,μ::Integer)=Derivative(d,μ)
Base.diff(d::Domain)=Base.diff(d,1)

Base.diff(d::Union(ProductDomain,TensorSpace),k::Integer)=Derivative(d,k)

# use conj(f.') for ArraySpace
Base.ctranspose(f::Fun)=differentiate(f)


∫(f::Fun)=integrate(f)
⨜(f::Fun)=cumsum(f)

for OP in (:Σ,:∮,:⨍,:⨎)
    @eval $OP(f::Fun)=sum(f)
end