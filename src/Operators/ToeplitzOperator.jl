export ToeplitzOperator, HankelOperator, LaurentOperator






type ToeplitzOperator{T<:Number} <: BandedOperator{T}
    negative::Vector{T}
    nonnegative::Vector{T}
end

ToeplitzOperator{T<:Number,Q<:Number}(V::Vector{T},W::Vector{Q})=ToeplitzOperator{promote_type(T,Q)}(V,W)
function ToeplitzOperator(V::Vector)
    W=V[2:end]
    V=copy(V)
    V[1]*=2
    ToeplitzOperator(W,V)
end





function toeplitz_addentries!(v::Vector,A,kr::Range)
    if !isempty(v)
        v1=v[1]
        M=maxabs(v)
        if abs(v1) > M*10eps(eltype(v))
            for k=kr
                A[k,k]+=2v1
            end
        end

        for j=2:length(v)
            vj=v[j]
            if abs(vj)> M*10eps(eltype(v))
                for k = kr
                    A[k,k+j-1]+=vj
                end
                for k = intersect(j:kr[end],kr)
                    A[k,k-j+1]+=vj
                end
            end
        end
    end
    A
end

function toeplitz_addentries!(v::Vector,A::BandedMatrix,kr::UnitRange)
    if !isempty(v)
        @inbounds v1=v[1]
        @simd for k=kr
            @inbounds A.data[A.l+1,k]+=2v1
        end

        for j=2:length(v)
            @inbounds vj=v[j]
            @simd for k = kr
                @inbounds A.data[j+A.l,k] +=vj
            end
            @simd for k = max(kr[1],j):kr[end]
                @inbounds A.data[2-j+A.l,k]+=vj
            end
        end
    end
    A
end

function toeplitz_addentries!(c::Number,neg::Vector,pos::Vector,A,kr::Range)
    for k=kr,j=1:min(length(neg),k-1)
        A[k,k-j] += c*neg[j]
    end
    for k=kr,j=1:length(pos)
        A[k,k+j-1] += c*pos[j]
    end

    A
end

toeplitz_addentries!(neg::Vector,pos::Vector,A,kr::Range)=toeplitz_addentries!(1,neg,pos,A,kr)



addentries!(T::ToeplitzOperator,A,kr::Range)=toeplitz_addentries!(1,T.negative,T.nonnegative,A,kr)
bandinds(T::ToeplitzOperator)=(-length(T.negative),length(T.nonnegative)-1)


## Hankel Operator


type HankelOperator{T<:Number} <: BandedOperator{T}
    coefficients::Vector{T}
end

HankelOperator(f::Fun)=HankelOperator(f.coefficients)

function hankel_addentries!(c::Number,v::Vector,A,kr::Range)
    M=maxabs(v)
    for j=1:length(v)
        vj=c*v[j]
        if abs(vj)>M*10eps(eltype(v))
            for k=intersect(kr,1:j)
                if j + 1 >= k+1
                    A[k,j-k+1] += vj
                end
            end
        end
    end

    A
end

hankel_addentries!(v::Vector,A,kr::Range)=hankel_addentries!(1,v,A,kr)


addentries!(T::HankelOperator,A,kr::Range)=hankel_addentries!(T.coefficients,A,kr)

bandinds(T::HankelOperator)=(1-length(T.coefficients),length(T.coefficients)-1)



## Laurent Operator

type LaurentOperator{T<:Number} <: BandedOperator{T}
    negative::Vector{T}
    nonnegative::Vector{T}
end



function shiftrowrange(kr::Range)
    if isodd(kr[1]) && isodd(kr[end])
        poskr=div(kr[1]-1,2):div(kr[end]-1,2)
        negkr=-div(kr[end]-1,2):-div(kr[1]+1,2)
    elseif isodd(kr[1]) # && iseven(kr[end])
        poskr=div(kr[1]-1,2):div(kr[end],2)-1
        negkr=-div(kr[end],2):-div(kr[1]+1,2)
    elseif isodd(kr[end]) # && iseven(kr[1])
        poskr=div(kr[1],2):div(kr[end]-1,2)
        negkr=-div(kr[end]-1,2):-div(kr[1],2)
    else # iseven(kr[end]) && iseven(kr[1])
        poskr=div(kr[1],2):div(kr[end],2)-1
        negkr=-div(kr[end],2):-div(kr[1],2)
    end
    negkr,poskr
end

# function laurent_addentries!(v::Vector,A,kr::Range)
#     negkr,poskr=shiftrowrange(kr)
#     br=1-length(v):length(v)-1
#
#     for k=poskr,j=br
#         if j≥0 # && k≥0
#             A[2k+1,2j+1] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         else # && k≥0
#             A[2k+1,-2j] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         end
#     end
#
#     for k=negkr,j=br
#         if j≥0 # && k <0
#             A[-2k,2j+1] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         else # j<0 && k<0
#             A[-2k,-2j] += (j ==0) ? 2v[1] : v[abs(j)+1]
#         end
#     end
#
#     A
# end




function laurent_addentries!(neg::Vector,pos::Vector,A,kr::Range)
    negkr,poskr=shiftrowrange(kr)

    for k=poskr
        for j=1:length(neg)
            if k-j≥0 # && k≥0
                #A[k,k+j]=v[j]
                # k->2k+1
                A[2k+1,2k-2j+1] += neg[j]
            else # k+j<0 && k≥0

                #A[k,k+j]=v[j]
                A[2k+1,-2k+2j] += neg[j]
            end
        end
        for j=0:length(pos)-1
            if k+j≥0 # && k≥0
                #A[k,k+j]=v[j]
                # k->2k+1
                A[2k+1,2k+2j+1] += pos[j+1]
            else # k+j<0 && k≥0

                #A[k,k+j]=v[j]
                A[2k+1,-2k-2j] += pos[j+1]
            end
        end
    end

    for k=negkr
        for j=1:length(neg)
            if k-j<0 # && k <0
                A[-2k,-2k+2j] += neg[j]
            else # k+j>0 && k<0
                A[-2k,2k-2j+1] += neg[j]
            end
        end
        for j=0:length(pos)-1
            if k+j<0 # && k <0
                A[-2k,-2k-2j] += pos[j+1]
            else # k+j>0 && k<0
                A[-2k,2k+2j+1] += pos[j+1]
            end
        end
    end

    A
end




addentries!(T::LaurentOperator,A,kr::Range)=laurent_addentries!(T.negative,T.nonnegative,A,kr)

shiftbandinds(T::LaurentOperator)=-length(T.negative),length(T.nonnegative)-1
function bandinds(T::LaurentOperator)
    sbi=shiftbandinds(T)
    min(2sbi[1],-2sbi[end]),max(2sbi[end],-2sbi[1])
end




## inv

function Base.inv(T::ToeplitzOperator)
    @assert length(T.nonnegative)==1
    ai=linsolve(T,[1.0];maxlength=100000)
    ToeplitzOperator(ai[2:end],ai[1:1])
end

function Fun(T::ToeplitzOperator)
   if length(T.nonnegative)==1
      Fun([T.nonnegative;T.negative],Taylor())
    elseif length(T.negative)==0
        Fun(T.nonnegative,Hardy{false}())
    else
        Fun(interlace(T.nonnegative,T.negative),Laurent(Circle()))
    end
end




