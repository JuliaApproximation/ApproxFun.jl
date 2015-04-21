
## boundary

∂{A<:Interval,B<:Interval}(d::ProductDomain{@compat(Tuple{A,B})})=PiecewiseInterval(d[1].a+im*d[2].a,d[1].b+im*d[2].a,d[1].b+im*d[2].b,d[1].a+im*d[2].b,d[1].a+im*d[2].a)
∂{A<:Interval,B<:PeriodicInterval}(d::ProductDomain{@compat(Tuple{A,B})})=UnionDomain([d[1].b+im*d[2],d[1].a+im*reverse(d[2])])
∂{A<:PeriodicInterval,B<:Interval}(d::ProductDomain{@compat(Tuple{A,B})})=UnionDomain([d[1]+im*d[2].a,reverse(d[1])+im*d[2].b])


## Union


