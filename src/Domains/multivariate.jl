
## boundary

∂{A<:Interval,B<:Interval}(d::ProductDomain{@compat(Tuple{A,B})})=PiecewiseInterval(d[1].a+im*d[2].a,d[1].b+im*d[2].a,d[1].b+im*d[2].b,d[1].a+im*d[2].b,d[1].a+im*d[2].a)
∂{A<:Interval,B<:PeriodicInterval}(d::ProductDomain{@compat(Tuple{A,B})})=UnionDomain([d[1].b+im*d[2],d[1].a+im*reverse(d[2])])
∂{A<:PeriodicInterval,B<:Interval}(d::ProductDomain{@compat(Tuple{A,B})})=UnionDomain([d[1]+im*d[2].a,reverse(d[1])+im*d[2].b])
∂{A<:PeriodicInterval,B<:PeriodicInterval}(d::ProductDomain{@compat(Tuple{A,B})})=EmptyDomain()


## Union



function Base.join{IT<:Interval}(p1::Vector{IT},p2::Vector{IT})
    for k=length(p1):-1:1,j=length(p2):-1:1
        if p1[k]==reverse(p2[j])
            deleteat!(p1,k)
            deleteat!(p2,j)
            break
        end
    end
    [p1;p2]
end

function ∂{PD<:ProductDomain}(d::UnionDomain{PD})
    bnd=map(d->pieces(∂(d)),d.domains)
    PiecewiseInterval(reduce(join,bnd))
end

