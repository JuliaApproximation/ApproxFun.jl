

Derivative(sp::PiecewiseSpace)=DerivativeWrapper(interlace(diagm(map(Derivative,sp.spaces))),0)
Derivative(sp::PiecewiseSpace,k::Integer)=DerivativeWrapper(interlace(diagm(map(s->Derivative(s,k),sp.spaces))),k)