immutable PiecewiseInterval{T<:Number} <:Domain{T}
    points::Vector{T}
end

PiecewiseInterval(d::Number...)=PiecewiseInterval([d...])
==(a::PiecewiseInterval,b::PiecewiseInterval)=a.points==b.points

Base.length(d::PiecewiseInterval)=length(d.points)-1
Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])
isperiodic(d::PiecewiseInterval)=first(d.points)==last(d.points)


function points(d::PiecewiseInterval,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:length(d)]...)]
end

pieces(d::PiecewiseInterval,k)=d[k]
pieces{T}(d::PiecewiseInterval{T})=Interval{T}[d[k] for k=1:length(d)]

Base.rand(d::PiecewiseInterval)=rand(d[rand(1:length(d))])
checkpoints{T}(d::PiecewiseInterval{T})=mapreduce(checkpoints,union,pieces(d))

for OP in (:(Base.first),:(Base.last))
    @eval $OP(d::PiecewiseInterval)=$OP(d.points)
end

Base.reverse(d::PeriodicInterval)=PeriodicInterval(d.b,d.a)
