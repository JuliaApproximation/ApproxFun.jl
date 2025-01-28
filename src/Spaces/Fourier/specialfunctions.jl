function Base.real{DD}(f::Fun{Laurent{DD}})
    n=ncoefficients(f)
    cfs=f.coefficients

    ret=Array(Float64,iseven(n)?n+1:n)
    ret[1]=real(cfs[1])

    for k=2:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]=imag(cfs[k])
        ret[k+1]=real(cfs[k])
    end
    for k=3:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]+=real(cfs[k])
        ret[k-1]-=imag(cfs[k])
    end


    Fun(ret,Fourier(domain(f)))
end


function Base.imag{DD}(f::Fun{Laurent{DD}})
    n=ncoefficients(f)
    cfs=f.coefficients

    ret=Array(Float64,iseven(n)?n+1:n)
    ret[1]=imag(cfs[1])

    for k=2:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]=-real(cfs[k])
        ret[k+1]=imag(cfs[k])
    end
    for k=3:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]+=imag(cfs[k])
        ret[k-1]+=real(cfs[k])
    end


    Fun(ret,Fourier(domain(f)))
end
