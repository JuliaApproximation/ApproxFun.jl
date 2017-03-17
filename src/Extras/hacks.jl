## Functions that depend on the structure of BandedMatrix


function pad!(A::BandedMatrix,n,m)
    A.data=pad(A.data,size(A.data,1),m)
    A.m=n
    A
end



# linear algebra


## Constructors that involve MultivariateFun

# function Fun(f::Function,d::BivariateSpace)
#     if f==zero
#         zeros(d)
#     elseif hasnumargs(f,2)
#         Fun(ProductFun(f,d))
#     else
#         Fun(ProductFun((x,y)->f((x,y)),d))
#     end
# end
# function Fun(f::Function,d::BivariateSpace,n::Integer)
#     if hasnumargs(f,2)
#         defaultFun(x->f(x...),d,n)
#     else
#         defaultFun(f,d,n)
#     end
# end

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



#TODO: Remove. This is a temporary fix while waiting for a pull request to be merged.
function Base.norm{N, T}(a::Vec{N, T}, p)
    isinf(p) && return maximum(abs,a)
    ret = abs(a[1])^p
    for k = 2:N
        ret += abs(a[k])^p
    end
    ret^(1/p)
end
