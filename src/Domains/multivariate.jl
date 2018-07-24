include("ProductDomain.jl")


## boundary

∂(d::ProductDomain{Tuple{A,B}}) where {A<:IntervalOrSegment,B<:IntervalOrSegment} =
    PiecewiseSegment([Vec(factor(d,1).a,factor(d,2).a),
                      Vec(factor(d,1).b,factor(d,2).a),
                      Vec(factor(d,1).b,factor(d,2).b),
                      Vec(factor(d,1).a,factor(d,2).b),
                      Vec(factor(d,1).a,factor(d,2).a)])
∂(d::ProductDomain{Tuple{A,B}}) where {A<:IntervalOrSegment,B<:PeriodicInterval} =
    UnionDomain((PeriodicInterval(Vec(factor(d,1).b,factor(d,2).a),Vec(factor(d,1).b,factor(d,2).b)),
        PeriodicInterval(Vec(factor(d,1).a,factor(d,2).b),Vec(factor(d,1).a,factor(d,2).a))))
∂(d::ProductDomain{Tuple{A,B}}) where {A<:PeriodicInterval,B<:IntervalOrSegment} =
    UnionDomain((PeriodicInterval(Vec(factor(d,1).a,factor(d,2).a),Vec(factor(d,1).b,factor(d,2).a)),
        PeriodicInterval(Vec(factor(d,1).b,factor(d,2).b),Vec(factor(d,1).a,factor(d,2).b))))
∂(d::ProductDomain{Tuple{A,B}}) where {A<:PeriodicInterval,B<:PeriodicInterval} = EmptyDomain()



## Union
function Base.join(p1::AbstractVector{IT},p2::AbstractVector{IT}) where IT<:IntervalOrSegment
    for k=length(p1):-1:1,j=length(p2):-1:1
        if p1[k]==reverse(p2[j])
            deleteat!(p1,k)
            deleteat!(p2,j)
            break
        end
    end
    [p1;p2]
end

function ∂(d::UnionDomain{Tuple{PD1,PD2}}) where {PD1<:ProductDomain,PD2<:ProductDomain}
    bnd=map(d->components(∂(d)),d.domains)
    PiecewiseSegment(reduce(join,bnd))
end
