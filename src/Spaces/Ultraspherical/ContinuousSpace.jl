

struct ContinuousSpace{T,R} <: Space{PiecewiseSegment{T},R}
    domain::PiecewiseSegment{T}
end

ContinuousSpace(d::PiecewiseSegment{T}) where {T} =
    ContinuousSpace{T,real(eltype(T))}(d)


Space(d::PiecewiseSegment) = ContinuousSpace(d)

isperiodic(C::ContinuousSpace) = isperiodic(domain(C))

spacescompatible(a::ContinuousSpace,b::ContinuousSpace) = domainscompatible(a,b)
conversion_rule(a::ContinuousSpace,
                b::PiecewiseSpace{CD,DD,RR}) where {CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},DD,RR<:Real} where {DDD,RRR} = a

plan_transform(sp::ContinuousSpace,vals::AbstractVector) =
    TransformPlan{eltype(vals),typeof(sp),false,Nothing}(sp,nothing)

function *(P::TransformPlan{T,SS,false},vals::AbstractVector{T}) where {T,SS<:ContinuousSpace}
    S = P.space
    n=length(vals)
    d=domain(S)
    K=ncomponents(d)
    k=div(n,K)

    PT=promote_type(prectype(d),eltype(vals))
    if k==0
        vals
    elseif isperiodic(d)
        ret=Array{PT}(undef, max(K,n-K))
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(component(d,j)),
                          vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
                ret[2]=cfs[1]+cfs[2]
            elseif j < K
                ret[j+1]=cfs[1]+cfs[2]
            end
            ret[K+j:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(component(d,j)),
                          vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
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
        ret=Array{PT}(undef, n-K+1)
        r=n-K*k

        for j=1:r
            cfs=transform(ChebyshevDirichlet{1,1}(component(d,j)),
                          vals[(j-1)*(k+1)+1:j*(k+1)])
            if j==1
                ret[1]=cfs[1]-cfs[2]
            end

            ret[j+1]=cfs[1]+cfs[2]
            ret[K+j+1:K:end]=cfs[3:end]
        end

        for j=r+1:K
            cfs=transform(ChebyshevDirichlet{1,1}(component(d,j)),
                          vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])

            if length(cfs) ≤ 1
                ret .= cfs
            else
                if j==1
                    ret[1]=cfs[1]-cfs[2]
                end
                ret[j+1]=cfs[1]+cfs[2]
                ret[K+j+1:K:end]=cfs[3:end]
            end
        end

        ret
    end
end

components(S::ContinuousSpace) = map(ChebyshevDirichlet{1,1},components(domain(S)))
canonicalspace(S::ContinuousSpace) = PiecewiseSpace(components(S))
convert(::Type{PiecewiseSpace}, S::ContinuousSpace) = canonicalspace(S)

blocklengths(C::ContinuousSpace) = Fill(ncomponents(C.domain),∞)

block(C::ContinuousSpace,k) = Block((k-1)÷ncomponents(C.domain)+1)


## components

components(f::Fun{CS},j::Integer) where {CS<:ContinuousSpace} = components(Fun(f,canonicalspace(f)),j)
components(f::Fun{CS}) where {CS<:ContinuousSpace} = components(Fun(f,canonicalspace(space(f))))


function points(f::Fun{CS}) where {CS<:ContinuousSpace}
    n=ncoefficients(f)
    d=domain(f)
    K=ncomponents(d)

    m=isperiodic(d) ? max(K,n+2K-1) : n+K
    points(f.space,m)
end

## Conversion

coefficients(cfsin::AbstractVector,A::ContinuousSpace,B::PiecewiseSpace) =
    defaultcoefficients(cfsin,A,B)

coefficients(cfsin::AbstractVector,A::ContinuousSpace,B::ContinuousSpace) =
    default_Fun(Fun(A,cfsin),B).coefficients


# We implemnt conversion between continuous space and PiecewiseSpace with Chebyshev dirichlet
function Conversion(ps::PiecewiseSpace{CD,DD,RR},cs::ContinuousSpace) where {CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},
                                                                    DD,RR<:Real} where {DDD,RRR}
    @assert ps == canonicalspace(cs)
    ConcreteConversion(ps,cs)
end

function Conversion(cs::ContinuousSpace,ps::PiecewiseSpace{CD,DD,RR}) where {CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},
                                                                    DD,RR<:Real} where {DDD,RRR}
    @assert ps == canonicalspace(cs)
    ConcreteConversion(cs,ps)
end


bandwidths(C::ConcreteConversion{PiecewiseSpace{CD,DD,RR},CS}) where {CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},
                                                                                    DD,RR<:Real,CS<:ContinuousSpace} where {DDD,RRR} =
    1,ncomponents(domain(rangespace(C)))


function getindex(C::ConcreteConversion{PiecewiseSpace{CD,DD,RR},CS,T},
                  k::Integer,j::Integer) where {T,DD,RR<:Real,CS<:ContinuousSpace,
                                                CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}}}  where {DDD,RRR}
    d=domain(rangespace(C))
    K=ncomponents(d)
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


bandwidths(C::ConcreteConversion{<:ContinuousSpace,
                                     PiecewiseSpace{CD,DD,RR}}) where {CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},
                                              DD,RR<:Real}  where {DDD,RRR} =
            isperiodic(domainspace(C)) ? (2ncomponents(domain(rangespace(C)))-1,1) :
                                         (ncomponents(domain(rangespace(C))),1)

function getindex(C::ConcreteConversion{<:ContinuousSpace,
                                            PiecewiseSpace{CD,DD,RR},T},
                      k::Integer,j::Integer) where {T,CD<:Tuple{Vararg{ChebyshevDirichlet{1,1,DDD,RRR}}},
                                        DD,RR<:Real} where {DDD,RRR}
    d=domain(domainspace(C))
    K=ncomponents(d)
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


Dirichlet(S::TensorSpace{<:Tuple{<:ChebyshevDirichlet{1,1,<:IntervalOrSegment{T},R1},
                                  <:ChebyshevDirichlet{1,1,<:IntervalOrSegment{T},R2}}},k) where {T,R1,R2} =
    k == 0 ? ConcreteDirichlet(S,0) : tensor_Dirichlet(S,k)

Dirichlet(d::ProductDomain{<:Tuple{<:IntervalOrSegment{T},<:IntervalOrSegment{T}}}) where {T<:Real} =
    Dirichlet(ChebyshevDirichlet{1,1}(factor(d,1))*ChebyshevDirichlet{1,1}(factor(d,2)))

isblockbanded(::Dirichlet{TensorSpace{Tuple{CD,CD},DD,RR}}) where {CD<:ChebyshevDirichlet,DD<:BivariateDomain,RR} =
    true

blockbandwidths(::Dirichlet{TensorSpace{Tuple{CD,CD},DD,RR}}) where {CD<:ChebyshevDirichlet,DD<:BivariateDomain,RR} =
    (0,2)

colstop(B::Dirichlet{TensorSpace{Tuple{CD,CD},DD,RR}},j::Integer) where {CD<:ChebyshevDirichlet,DD<:BivariateDomain,RR} =
    j ≤ 3 ? 4 : 4(block(domainspace(B),j).n[1]-1)


function getindex(B::ConcreteDirichlet{TensorSpace{Tuple{CD,CD},DD,RR}},
k::Integer,j::Integer) where {CD<:ChebyshevDirichlet,DD<:BivariateDomain,RR}
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
        K = Int(block(rs,k))
        J = Int(block(ds,j))
        m = mod(k-1,4)
        s,t =  blockstart(ds,J),  blockstop(ds,J)
        if K == J-1 && (m == 1  && j == s ||
                       (m == 0  && j == t))
            one(T)
        elseif K == J-1 && ((m == 3 && j == s) ||
                            (m == 2 && j == t))
            iseven(K) ? one(T) : -one(T)
        elseif K == J-2 && m == 1 && j == s+1
            one(T)
        elseif K == J-2 && m == 2 && j == t-1
            iseven(K) ? one(T) : -one(T)
        elseif K == J-2 && m == 0 && j == t-1
            -one(T)
        elseif K == J-2 && m == 3 && j == s+1
            iseven(K) ? -one(T) : one(T)
        else
            zero(T)
        end
    end
end


function BlockBandedMatrix(S::SubOperator{T,ConcreteDirichlet{TensorSpace{Tuple{CD,CD},DD,RR},
                                                    CSP,TT},
                                Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,CD<:ChebyshevDirichlet,DD<:BivariateDomain,RR,CSP,TT}
    P=parent(S)
    ret=BlockBandedMatrix(Zeros, S)
    kr,jr=parentindices(S)

    K1=block(rangespace(P),kr[1])
    Kr1=blockstart(rangespace(P),K1)
    J1=block(domainspace(P),jr[1])
    Jr1=blockstart(domainspace(P),J1)

    if ret.rows[1]>0 && ret.cols[1]==1
        view(ret,Block(1),Block(1))[:,1] = 1
    end


    if ret.rows[1] > 0 && ret.cols[2] > 0
        B=view(ret,Block(1),Block(2))

        k_sh = kr[1]-1; j_sh = max(jr[1]-2,0)
        if j_sh == 0
            # first column
            k_sh == 0 && (B[1,1-j_sh]=-1)
            k_sh ≤  1 && (B[2-k_sh,1-j_sh]=-1)
            k_sh ≤ 2 && (B[3-k_sh,1-j_sh]=1)
            B[4-k_sh,1-j_sh]=1
        end
        # second column
        k_sh == 0 && (B[1-k_sh,2-j_sh]=-1)
        k_sh ≤ 1  && (B[2-k_sh,2-j_sh]=1)
        k_sh ≤ 2  && (B[3-k_sh,2-j_sh]=1)

        B[4-k_sh,2-j_sh]=-1
    end


    if ret.rows[1] > 0 && ret.cols[3] > 0
        B=view(ret,Block(1),Block(3))

        k_sh = kr[1]-1; j_sh = max(jr[1]-4,0)
        # second column
        k_sh == 0 && (B[1-k_sh,2-j_sh]=1)
        k_sh ≤ 1 && (B[2-k_sh,2-j_sh]=-1)
        k_sh ≤ 2 && (B[3-k_sh,2-j_sh]=1)
        B[4-k_sh,2-j_sh]=-1
    end
    for K=Block(2):2:Block(min(length(ret.rows),length(ret.cols)-1))
        J = K+1  # super-diagonal block
        N = ret.rows[K.n[1]]
        M = ret.cols[J.n[1]]
        if N ≠ 0 && M ≠ 0
            # calculate shift
            k_sh = K == K1 ? kr[1]-Kr1 : 0
            j_sh = J == J1 ? jr[1]-Jr1 : 0
            B = view(ret,K,J)

            1 ≤ 2-k_sh ≤ N && j_sh == 0 && (B[2-k_sh,1-j_sh]=1)
            1 ≤ 4-k_sh ≤ N && j_sh == 0 && (B[4-k_sh,1-j_sh]=1)
            k_sh == 0 && 1 ≤ J.n[1]-j_sh ≤ M && (B[1-k_sh,J.n[1]-j_sh]=1)
            k_sh ≤ 2 &&  1 ≤ J.n[1]-j_sh ≤ M && (B[3-k_sh,J.n[1]-j_sh]=1)
        end
    end
    for K=Block(3):2:Block(min(length(ret.rows),length(ret.cols)-1))
        J = K+1  # super-diagonal block
        N = ret.rows[K.n[1]]
        M = ret.cols[J.n[1]]
        if N ≠ 0 && M ≠ 0
            # calculate shift
            k_sh = K == K1 ? kr[1]-Kr1 : 0
            j_sh = J == J1 ? jr[1]-Jr1 : 0
            B = view(ret,K,J)

            1 ≤ 2-k_sh ≤ N && j_sh == 0 && (B[2-k_sh,1-j_sh]=1)
            1 ≤ 4-k_sh ≤ N && j_sh == 0 && (B[4-k_sh,1-j_sh]=-1)
            k_sh == 0 && 1 ≤ J.n[1]-j_sh ≤ M && (B[1-k_sh,J.n[1]-j_sh]=1)
            1 ≤ 3-k_sh ≤ N &&  1 ≤ J.n[1]-j_sh ≤ M && (B[3-k_sh,J.n[1]-j_sh]=-1)
        end
    end
    for K=Block(2):2:Block(min(length(ret.rows),length(ret.cols)-2))
        J = K+2  # super-diagonal block
        N = ret.rows[K.n[1]]
        M = ret.cols[J.n[1]]

        if N ≠ 0 && M ≠ 0
            B=view(ret,K,J)
            # calculate shift
            k_sh = K == K1 ? kr[1]-Kr1 : 0
            j_sh = J == J1 ? jr[1]-Jr1 : 0
            B = view(ret,K,J)

            1 ≤ 2-k_sh ≤ N && 1 ≤ 2-j_sh ≤ M && (B[2-k_sh,2-j_sh]=1)
            1 ≤ 4-k_sh ≤ N && 1 ≤ 2-j_sh ≤ M && (B[4-k_sh,2-j_sh]=-1)
            k_sh == 0 && 1 ≤ J.n[1]-j_sh-1 ≤ M && (B[1,J.n[1]-j_sh-1]=-1)
            1 ≤ 3-k_sh ≤ N &&  1 ≤ J.n[1]-j_sh-1 ≤ M && (B[3-k_sh,J.n[1]-j_sh-1]=1)
        end
    end
    for K=Block(3):2:Block(min(length(ret.rows),length(ret.cols)-2))
        J = K+2
        B=view(ret,K,J)
        N = ret.rows[K.n[1]]
        M = ret.cols[J.n[1]]

        if N ≠ 0 && M ≠ 0
            B=view(ret,K,J)
            # calculate shift
            k_sh = K == K1 ? kr[1]-Kr1 : 0
            j_sh = J == J1 ? jr[1]-Jr1 : 0
            B = view(ret,K,J)

            1 ≤ 2-k_sh ≤ N && 1 ≤ 2-j_sh ≤ M && (B[2-k_sh,2-j_sh]=1)
            1 ≤ 4-k_sh ≤ N && 1 ≤ 2-j_sh ≤ M && (B[4-k_sh,2-j_sh]=1)
            k_sh == 0 && 1 ≤ J.n[1]-j_sh-1 ≤ M && (B[1,J.n[1]-j_sh-1]=-1)
            1 ≤ 3-k_sh ≤ N &&  1 ≤ J.n[1]-j_sh-1 ≤ M && (B[3-k_sh,J.n[1]-j_sh-1]=-1)
        end
    end

    ret
end


union_rule(A::PiecewiseSpace, B::ContinuousSpace) = union(A, convert(PiecewiseSpace, B))
union_rule(A::ConstantSpace, B::ContinuousSpace) = B

function approx_union(a::AbstractVector{T}, b::AbstractVector{V}) where {T,V}
    ret = sort!(union(a,b))
    for k in length(ret)-1:-1:1
        isapprox(ret[k] , ret[k+1]; atol=10eps()) && deleteat!(ret, k+1)
    end
    ret
end



function union_rule(A::ContinuousSpace{<:Real}, B::ContinuousSpace{<:Real})
    p_A,p_B = domain(A).points, domain(B).points
    a,b = minimum(p_A),  maximum(p_A)
    c,d = minimum(p_B),  maximum(p_B)
    @assert !isempty((a..b) ∩ (c..d))
    ContinuousSpace(PiecewiseSegment(approx_union(p_A, p_B)))
end
