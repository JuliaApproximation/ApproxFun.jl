

for TYP in (:ReverseOrientation,:Reverse)
    WRAP = Meta.parse(string(TYP)*"Wrapper")
    @eval begin
        abstract type $TYP{T} <: Operator{T} end

        struct $WRAP{OS,T} <: Operator{T}
            op::OS
        end

        $WRAP(op::Operator) = $WRAP{typeof(op),eltype(op)}(op)
        convert(::Type{Operator{T}},op::$TYP) where {T} = $TYP{T}()
        convert(::Type{Operator{T}},op::$WRAP) where {T} = $WRAP(Operator{T}(op.op))::Operator{T}

        @wrapper $WRAP
    end
end
