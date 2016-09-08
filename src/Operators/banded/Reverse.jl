

for TYP in (:ReverseOrientation,:Reverse)
    WRAP = parse(string(TYP)*"Wrapper")
    @eval begin
        abstract $TYP{T} <: Operator{T}

        immutable $WRAP{OS,T} <: Operator{T}
            op::OS
        end

        @wrapper $WRAP
    end
end
