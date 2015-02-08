

function spaceconversion(f::Vector,sp1::JacobiWeight{Chebyshev},sp2::JacobiWeight{Chebyshev})
    α,β=sp1.α,sp1.β
    c,d=sp2.α,sp2.β
    if c==α && d==β
        f
    elseif c>α && d>β
        spaceconversion(divide_singularity(f),JacobiWeight(α+1,β+1,sp1.space),sp2)
    elseif c>α
        spaceconversion(divide_singularity(-1,f),JacobiWeight(α+1,β,sp1.space),sp2)
    elseif d>β
        spaceconversion(divide_singularity(1,f),JacobiWeight(α,β+1,sp1.space),sp2)
    else
        error("Need to implement decreasing jacobi")
    end
end



transform(sp::JacobiWeight{Chebyshev},vals::Vector)=chebyshevtransform(vals./jacobiweight(sp,points(sp,length(vals)));kind=1)
itransform(sp::JacobiWeight{Chebyshev},cfs::Vector)=ichebyshevtransform(cfs).*jacobiweight(sp,points(sp,length(cfs));kind=1)



