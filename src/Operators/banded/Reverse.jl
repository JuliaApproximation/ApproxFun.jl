

for TYP in (:ReverseOrientation,:Reverse)
    WRAP = parse(string(TYP)*"Wrapper")
    @eval begin
        abstract type $TYP{T} <: Operator{T} end

        struct $WRAP{OS,T} <: Operator{T}
            op::OS
        end

        $WRAP(op::Operator) = $WRAP{typeof(op),eltype(op)}(op)
        convert{T}(::Type{Operator{T}},op::$TYP) = $TYP{T}()
        convert{T}(::Type{Operator{T}},op::$WRAP) = $WRAP(Operator{T}(op.op))::Operator{T}

        @wrapper $WRAP
    end
end
