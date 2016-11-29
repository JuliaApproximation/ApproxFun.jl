include("ProductDomain.jl")


## boundary

∂{A<:Segment,B<:Segment}(d::ProductDomain{Tuple{A,B}}) =
    PiecewiseInterval([Vec(d[1].a,d[2].a),
                      Vec(d[1].b,d[2].a),
                      Vec(d[1].b,d[2].b),
                      Vec(d[1].a,d[2].b),
                      Vec(d[1].a,d[2].a)])
∂{A<:Segment,B<:PeriodicInterval}(d::ProductDomain{Tuple{A,B}}) =
    UnionDomain((PeriodicInterval(Vec(d[1].b,d[2].a),Vec(d[1].b,d[2].b)),
        PeriodicInterval(Vec(d[1].a,d[2].b),Vec(d[1].a,d[2].a))))
∂{A<:PeriodicInterval,B<:Segment}(d::ProductDomain{Tuple{A,B}}) =
    UnionDomain((PeriodicInterval(Vec(d[1].a,d[2].a),Vec(d[1].b,d[2].a)),
        PeriodicInterval(Vec(d[1].b,d[2].b),Vec(d[1].a,d[2].b))))
∂{A<:PeriodicInterval,B<:PeriodicInterval}(d::ProductDomain{Tuple{A,B}}) = EmptyDomain()



## Union



function Base.join{IT<:Segment}(p1::AbstractVector{IT},p2::AbstractVector{IT})
    for k=length(p1):-1:1,j=length(p2):-1:1
        if p1[k]==reverse(p2[j])
            deleteat!(p1,k)
            deleteat!(p2,j)
            break
        end
    end
    [p1;p2]
end

function ∂{PD1<:ProductDomain,PD2<:ProductDomain}(d::UnionDomain{Tuple{PD1,PD2}})
    bnd=map(d->pieces(∂(d)),d.domains)
    PiecewiseInterval(reduce(join,bnd))
end
