

for TYP in (:ReverseOrientation,:Reverse)
    WRAP = parse(string(TYP)*"Wrapper")
    @eval begin
        abstract $TYP{T} <: Operator{T}

        immutable $WRAP{OS,T} <: Operator{T}
            op::OS
        end

        $WRAP(op::Operator) = $WRAP{typeof(op),eltype(op)}(op)

        @wrapper $WRAP
    end
end
