immutable Point{T} <: Domain{T,0}
    x::T
end



for op in (:*,:+,:-,:.*,:.+,:.-,:.^)
    @eval begin
        $op(c::Number,d::Point)=Point($op(c,d.x))
        $op(d::Point,c::Number)=Point($op(d.x,c))
    end
end

for op in (:/,:./)
    @eval $op(d::Point,c::Number)=Point($op(d.x,c))
end
