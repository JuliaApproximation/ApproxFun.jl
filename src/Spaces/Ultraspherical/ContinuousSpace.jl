

immutable ContinuousSpace{T} <: Space{RealBasis,PiecewiseInterval{T},1}
    domain::PiecewiseInterval{T}
end



Space(d::PiecewiseInterval) = ContinuousSpace(d)

isperiodic(C::ContinuousSpace) = isperiodic(domain(C))

spacescompatible(a::ContinuousSpace,b::ContinuousSpace) = domainscompatible(a,b)
conversion_rule{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}}}(a::ContinuousSpace,
                                                            b::PiecewiseSpace{CD,RealBasis}) = a

function transform(S::ContinuousSpace,vals::Vector)
    n=length(vals)
    d=domain(S)
    K=numpieces(d)
    k=div(n,K)

    PT=promote_type(eltype(eltype(d)),eltype(vals))
    if k==0
        vals
    elseif isperiodic(d)
        ret=Array(PT,max(K,n-K))
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
                ret[2]=cfs[1]+cfs[2]
            elseif j < K
                ret[j+1]=cfs[1]+cfs[2]
            end
            ret[K+j:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            if length(cfs)==1 && j <K
                ret[j+1]=cfs[1]
            elseif j==1
                ret[1]=cfs[1]-cfs[2]
                ret[2]=cfs[1]+cfs[2]
            elseif j < K
                ret[j+1]=cfs[1]+cfs[2]
            end
            ret[K+j:K:end]=cfs[3:end]
        end

        ret
    else
        ret=Array(PT,n-K+1)
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
            end

            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            if j==1
                ret[1]=cfs[1]-cfs[2]
            end
            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        ret
    end
end

canonicalspace(S::ContinuousSpace) = PiecewiseSpace(map(ChebyshevDirichlet{1,1},pieces(domain(S))))





blocklengths(C::ContinuousSpace) = repeated(numpieces(C.domain))

block(C::ContinuousSpace,k) = (k-1)÷numpieces(C.domain)+1


## pieces

Base.vec{T}(f::Fun{ContinuousSpace{T}},j::Integer) = vec(Fun(f,canonicalspace(f)),j)
Base.vec{T}(f::Fun{ContinuousSpace{T}}) = vec(Fun(f,canonicalspace(space(f))))
pieces{T}(f::Fun{ContinuousSpace{T}}) = vec(f)


function points{T}(f::Fun{ContinuousSpace{T}})
    n=ncoefficients(f)
    d=domain(f)
    K=numpieces(d)

    m=isperiodic(d)?max(K,n+2K-1):n+K
    points(f.space,m)
end

## Conversion

coefficients(cfsin::Vector,A::ContinuousSpace,B::PiecewiseSpace) =
    defaultcoefficients(cfsin,A,B)


# We implemnt conversion between continuous space and PiecewiseSpace with Chebyshev dirichlet
Conversion{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},
           DD}(ps::PiecewiseSpace{CD,RealBasis,DD,1},cs::ContinuousSpace) =
                ConcreteConversion(ps,cs)

Conversion{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},
           DD}(cs::ContinuousSpace,ps::PiecewiseSpace{CD,RealBasis,DD,1}) =
                ConcreteConversion(cs,ps)


bandinds{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},
         DD,T}(C::ConcreteConversion{PiecewiseSpace{CD,RealBasis,DD,1},ContinuousSpace{T}}) =
    -1,numpieces(domain(rangespace(C)))


function getindex{T,DD,TT,
                  CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}}}(C::ConcreteConversion{PiecewiseSpace{CD,RealBasis,DD,1},ContinuousSpace{TT},T},
                                                              k::Integer,j::Integer)
    d=domain(rangespace(C))
    K=numpieces(d)
    if isperiodic(d)
        if k==j==1
            one(T)
        elseif k==1 && j==K+1
            -one(T)
        elseif 2≤k≤K && (j==k-1 || j==K+k-1)
            one(T)
        elseif K<k && j==k+K
            one(T)
        else
            zero(T)
        end
    else
        if k==j==1
            one(T)
        elseif k==1 && j==K+1
            -one(T)
        elseif 2≤k≤K+1 && (j==k-1 || j==K+k-1)
            one(T)
        elseif K+1<k && j==k+K-1
            one(T)
        else
            zero(T)
        end
    end
end


bandinds{CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},
         DD,T}(C::ConcreteConversion{ContinuousSpace{T},
                                     PiecewiseSpace{CD,RealBasis,DD,1}}) =
            isperiodic(domainspace(C)) ? (1-2numpieces(domain(rangespace(C))),1) :
                                         (-numpieces(domain(rangespace(C))),1)

function getindex{T,CD<:Tuple{Vararg{ChebyshevDirichlet{1,1}}},TT,
                  DD}(C::ConcreteConversion{ContinuousSpace{TT},
                                            PiecewiseSpace{CD,RealBasis,DD,1},T},
                      k::Integer,j::Integer)
    d=domain(domainspace(C))
    K=numpieces(d)
    if isperiodic(d)
        if k < K && (j==k || j==k+1)
            one(T)/2
        elseif k==K && (j==K || j==1)
            one(T)/2
        elseif K+1≤k≤2K && j==k-K
            -one(T)/2
        elseif K+1≤k<2K && j==k-K+1
            one(T)/2
        elseif k==2K && j==1
            one(T)/2
        elseif k>2K && j==k-K
            one(T)
        else
            zero(T)
        end
    else
        if k≤K && (j==k || j==k+1)
            one(T)/2
        elseif K+1≤k≤2K && j==k-K
            -one(T)/2
        elseif K+1≤k≤2K && j==k-K+1
            one(T)/2
        elseif k>2K && j==k-K+1
            one(T)
        else
            zero(T)
        end
    end
end



# Dirichlet for Squares


Dirichlet{T}(S::TensorSpace{ChebyshevDirichlet{1,1,Interval{T}},ChebyshevDirichlet{1,1,Interval{T}}}) =
    ConcreteDirichlet(S)

Dirichlet{T<:Real}(d::ProductDomain{Tuple{Interval{T},Interval{T}}}) =
    Dirichlet(ChebyshevDirichlet{1,1}(d[1])*ChebyshevDirichlet{1,1}(d[2]))

isblockbanded{CD<:ChebyshevDirichlet,RB}(::Dirichlet{TensorSpace{Tuple{CD,CD},RB,2}}) =
    true

blockbandinds{CD<:ChebyshevDirichlet,RB}(::Dirichlet{TensorSpace{Tuple{CD,CD},RB,2}}) =
    (0,2)

colstop{CD<:ChebyshevDirichlet,RB}(B::Dirichlet{TensorSpace{Tuple{CD,CD},RB,2}},j::Integer) =
    j ≤ 3 ? 4 : 4(block(domainspace(B),j)-1)


function getindex{CD<:ChebyshevDirichlet,RB}(B::Dirichlet{TensorSpace{Tuple{CD,CD},RB,2}},
                                             k::Integer,j::Integer)
    T = eltype(B)
    ds = domainspace(B)
    rs = rangespace(B)
    if j == 1 && k ≤ 4
        one(T)
    elseif j == 2 && k ≤ 2
        -one(T)
    elseif j == 2 && k ≤ 4
        one(T)
    elseif j == 3 && (k == 1 || k == 4)
        -one(T)
    elseif j == 3 && (k == 2 || k == 3)
        one(T)
    elseif j == 5 && (k == 2 || k == 4)
        -one(T)
    elseif j == 5 && (k == 1 || k == 3)
        one(T)
    elseif j == 5 || j ≤ 3
        zero(T)
    else
        K = block(rs,k)
        J = block(ds,j)
        m = mod(k-1,4)
        s,t =  blockstart(ds,J),  blockstop(ds,J)
        if K == J-1 && (m == 1  && j == s ||
                       (m == 0  && j == t))
            one(T)
        elseif K == J-1 && ((m == 3 && j == s) ||
                            (m == 2 && j == t))
            iseven(K)?one(T):-one(T)
        elseif K == J-2 && m == 1 && j == s+1
            one(T)
        elseif K == J-2 && m == 2 && j == t-1
            iseven(K)?one(T):-one(T)
        elseif K == J-2 && m == 0 && j == t-1
            -one(T)
        elseif K == J-2 && m == 3 && j == s+1
            iseven(K)?-one(T):one(T)
        else
            zero(T)
        end
    end
end
