## Functions that depend on the structure of BandedMatrix


function pad!(A::BandedMatrix,n,m)
    A.data = pad(A.data,size(A.data,1),m)
    A.raxis = Base.OneTo(n)
    A
end



# linear algebra


## Constructors that involve MultivariateFun
Fun(f::Fun) = f # Fun of Fun should be like a conversion

function Fun(fin::Function)
    f = dynamic(fin)

    if hasnumargs(f,1)
        # check for tuple
        try
            f(0)
        catch ex
            if isa(ex,BoundsError)
                # assume its a tuple
                return Fun(f,ChebyshevInterval()^2)
            else
                throw(ex)
            end
        end

        Fun(f,ChebyshevInterval())
    elseif hasnumargs(f,2)
            Fun(f,ChebyshevInterval()^2)
    else
        error("Function not defined on interval or square")
    end
end




## These hacks support PDEs with block matrices



## dot for AbstractVector{Number} * AbstractVector{Fun}

function dot(c::AbstractVector{T},f::AbstractVector{F}) where {T<:Union{Number,Fun,MultivariateFun},F<:Union{Fun,MultivariateFun}}
    @assert length(c)==length(f)
    ret = conj(first(c))*first(f)
    for k = 2:length(c)
        ret += conj(c[k])*f[k]
    end
    ret
end


function dotu(c::AbstractVector{T},f::AbstractVector{F}) where {T<:Union{Fun,MultivariateFun,Number},F<:Union{Fun,MultivariateFun,Number}}
    @assert length(c)==length(f)
    isempty(c) && return zero(Base.promote_op(*,T,F))
    ret = c[1]*f[1]
    for k = 2:length(c)
        ret += c[k]*f[k]
    end
    ret
end



## Gets blockbandwidths working for SpectralMeasures

# TODO: Is this a good definition?
blockbandwidths(T::FiniteOperator{AT,TT,SS1,SS2}) where {AT,TT,SS1<:Union{EuclideanSpace,SequenceSpace},
              SS2<:Union{EuclideanSpace,SequenceSpace}} = bandwidths(T)
