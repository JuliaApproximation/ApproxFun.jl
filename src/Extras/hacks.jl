## Functions that depend on the structure of BandedMatrix


function pad!(A::BandedMatrix,n,m)
    A.data=pad(A.data,size(A.data,1),m)
    A.m=n
    A
end



# linear algebra


# function toeplitz_addentries!(v::Vector,A::BandedMatrix,kr::UnitRange)
#     if !isempty(v)
#         @inbounds v1=v[1]
#         @simd for k=kr
#             @inbounds A.data[A.u+1,k]+=2v1
#         end
#
#         for j=2:length(v)
#             @inbounds vj=v[j]
#             @simd for k = kr
#                 @inbounds A.data[A.u-j+2,k+j-1] +=vj
#             end
#             @simd for k = max(kr[1],j):kr[end]
#                 @inbounds A.data[A.u+j,k-j+1]+=vj
#             end
#         end
#     end
#     A
# end

## Constructors that involve MultivariateFun

function Fun(f::Function,d::BivariateSpace)
    if f==zero
        zeros(d)
    elseif hasnumargs(f,2)
        Fun(ProductFun(f,d))
    else
        Fun(ProductFun((x,y)->f((x,y)),d))
    end
end
function Fun(f::Function,d::BivariateSpace,n::Integer)
    if hasnumargs(f,2)
        defaultFun(x->f(x...),d,n)
    else
        defaultFun(f,d,n)
    end
end

function Fun(f::Function)
    if hasnumargs(f,1)
        # check for tuple
        try
            f(0)
        catch ex
            if isa(ex,BoundsError)
                # assume its a tuple
                return Fun(f,Interval()^2)
            else
                throw(ex)
            end
        end

        Fun(f,Interval())
    elseif hasnumargs(f,2)
            Fun(f,Interval()^2)
    else
        error("Function not defined on interval or square")
    end
end




## These hacks support PDEs with block matrices



## dot for vector{Number} * Vector{Fun}

function Base.dot{T<:Union{Number,Fun,MultivariateFun},F<:Union{Fun,MultivariateFun}}(c::Vector{T},f::Vector{F})
    @assert length(c)==length(f)
    ret=conj(first(c))*first(f)
    for k=2:length(c)
        ret+=conj(c[k])*f[k]
    end
    ret
end


# for TYP in (:Real,:Number)
#     @eval begin
#         function dotu{T<:$TYP,F<:Union{Fun,MultivariateFun}}(c::Vector{T},f::Vector{F})
#             @assert length(c)==length(f)
#             ret=first(c)*first(f)
#             for k=2:length(c)
#                 ret+=c[k]*f[k]
#             end
#             ret
#         end
#     end
# end

function dotu{T<:Union{Fun,MultivariateFun,Number},F<:Union{Fun,MultivariateFun,Number}}(c::Vector{T},f::Vector{F})
    @assert length(c)==length(f)
    ret=first(c)*first(f)
    for k=2:length(c)
        ret+=c[k]*f[k]
    end
    ret
end
