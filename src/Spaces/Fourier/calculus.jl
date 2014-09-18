##Differentiation and integration


differentiate(f::Fun{LaurentSpace})=Fun(interlace(fourierdiff(domain(f),deinterlace(f.coefficients))),f.space)
Base.sum(f::Fun{LaurentSpace})=fouriersum(domain(f),deinterlace(f.coefficients))
integrate(f::Fun{LaurentSpace})=Fun(interlace(fourierintegrate(domain(f),deinterlace(f.coefficients))),f.space)


fourierdiff(d::PeriodicInterval,cfs::ShiftVector)=tocanonicalD(d,0)*ShiftVector(1.im*[firstindex(cfs):-1],1.im*[0:lastindex(cfs)]).*cfs




function fourierintegrate(d::PeriodicInterval,cfs::ShiftVector)
    tol = 10eps()
    @assert abs(cfs[0]) < tol
    
    ##TODO: mapped domains
    
    @assert d.a ==-π
    @assert d.b ==π        
    ShiftVector(-1.im./[firstindex(cfs):-1],
                [0,(-1.im./[1:lastindex(cfs)])])
end

fouriersum(d::PeriodicInterval,cfs::ShiftVector)=cfs[0].*length(d)



function fourierdiff(d::Circle,cfs::ShiftVector)
        ##TODO: general radii
        @assert d.radius == 1.
        @assert d.center == 0

        # Now shift everything by one
        ShiftVector(
                        [([firstindex(cfs):-1].*cfs[firstindex(cfs):-1]),0],
                        [1:lastindex(cfs)].*cfs[1:lastindex(cfs)]
                        )
end



function fourierintegrate(d::Circle,cfs::ShiftVector)
    tol = 10eps()
    @assert abs(cfs[-1]) < tol        
    ##TODO: general radii        
    @assert d.radius == 1.
    @assert d.center == 0        
    
    # Now shift everything by one
    ShiftVector(
                    [cfs[firstindex(cfs):-1]./[firstindex(cfs):-1]],
                    [0,(cfs[0:lastindex(cfs)]./[1:lastindex(cfs)+1])]
                    )
end


function fouriersum{T}(d::Circle,cfs::ShiftVector{T})
    @assert d.radius == 1.
    @assert d.center == 0   
    if firstindex(cfs) <= -1
        cfs[-1]
    else
        zero(T)
    end
end

