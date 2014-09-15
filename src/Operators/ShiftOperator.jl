## BandedShiftOperator allows automatic interlacing

abstract BandedShiftOperator{T} <: BandedOperator{T}
abstract ShiftFunctional{T} <: Functional{T}

bandrange(b::BandedShiftOperator)=Range1(bandinds(b)...)