using ApproxFun
    import ApproxFun:BandedOperator,bandinds,DiagonalOperator,addentries!,domainspace,rangespace

immutable Reflection{S,T} <: DiagonalOperator{T}
    space::S
end

Reflection(S::FunctionSpace)=Reflection{typeof(S),Float64}(S)
domainspace(R::Reflection)=R.space
rangespace(R::Reflection)=R.space

function addentries!(R::Reflection,A,kr::Range)
   for k in kr
       A[k,k]+=(iseven(k)?-1:1)
   end
   A
end

f=Fun(exp)
R=Reflection(space(f))

(R*f)[.1]
f[-.1]

R[1:10,1:10]
D[1:10,1:10]

