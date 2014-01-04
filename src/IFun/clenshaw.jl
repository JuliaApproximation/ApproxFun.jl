function clenshaw(c,x)
    if maximum(abs(imag(x))) < 10*eps()
        clenshaw(c,real(x))
    else
        clenshaw(c,real(x))+im*clenshaw(c,imag(x))    
    end
end

function clenshaw(c::Vector,x::Real)

    x = 2x;
    bk1 = 0.0;
    bk2 = 0.0;
    for k = length(c):-1:2
        bk2, bk1 = bk1, c[k] + x * bk1 - bk2
    end

    c[1] + 0.5 * x * bk1 - bk2
end



function clenshaw{T<:Number,M<:Real}(c::Vector{T},x::Vector{M})
    @assert length(c) > 0 
    
    n = length(x)
    clenshaw(c,x,Array(T,n),Array(T,n),Array(T,n))
end


#Clenshaw routine for many Funs, x is a vector of same number of funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Vector{T})=clenshaw(c,x,Array(T,size(c)[2]),Array(T,size(c)[2]),Array(T,size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Vector{T},bk::Vector{T},bk1::Vector{T},bk2::Vector{T})
    n=size(c)[2] #number of funs
    m=size(c)[1] #number of coefficients

    bk1_v = unsafe_view(bk1)
    bk2_v = unsafe_view(bk2)
    bk_v = unsafe_view(bk)
    c_v = unsafe_view(c)
    x_v = unsafe_view(x)
    
    for i = 1:n
        bk1_v[i] = 0.
        bk2_v[i] = 0.
    end

    for k=m:-1:2
        for j=1:n
            ck = c_v[k,j]

            bk_v[j] = ck + 2x[j] * bk1_v[j] - bk2_v[j]
        end
        
        bk2_v, bk1_v, bk_v = bk1_v, bk_v, bk2_v
    end


    for j = 1:n
        ce = c_v[1,j]
        bk[j] = ce + x[j] * bk1_v[j] - bk2_v[j]
    end
    
    bk    
end



#Clenshaw routine for many Funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Real)=clenshaw(c,x,Array(T,size(c)[2]),Array(T,size(c)[2]),Array(T,size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Real,bk::Vector{T},bk1::Vector{T},bk2::Vector{T})
    n=size(c)[2] #number of funs
    m=size(c)[1] #number of coefficients

    bk1_v = unsafe_view(bk1)
    bk2_v = unsafe_view(bk2)
    bk_v = unsafe_view(bk)
    c_v = unsafe_view(c)
    
    for i = 1:n
        bk1_v[i] = 0.
        bk2_v[i] = 0.
    end

    for k=m:-1:2
        for j=1:n
            ck = c_v[k,j]

            bk_v[j] = ck + 2x * bk1_v[j] - bk2_v[j]
        end
        
        bk2_v, bk1_v, bk_v = bk1_v, bk_v, bk2_v
    end


    for j = 1:n
        ce = c_v[1,j]
        bk[j] = ce + x * bk1_v[j] - bk2_v[j]
    end
    
    bk    
end


#Clenshaw routine for many points
#Note that bk1, bk2, and bk are overwritten
function clenshaw{T<:Number,M<:Real}(c::Vector{T},x::Vector{M},bk::Vector{T},bk1::Vector{T},bk2::Vector{T})
    n = length(x)
#    x=2x

    x_v = unsafe_view(x)
    bk1_v = unsafe_view(bk1)
    bk2_v = unsafe_view(bk2)
    bk_v = unsafe_view(bk)
    c_v = unsafe_view(c)
    
    for i = 1:n
        bk1_v[i] = 0.
        bk2_v[i] = 0.
    end

    for k in  length(c):-1:2
        ck = c_v[k]
        for i in 1 : n
            bk_v[i] = ck + 2x_v[i] * bk1_v[i] - bk2_v[i]
        end
        bk2_v, bk1_v, bk_v = bk1_v, bk_v, bk2_v
    end

    ce = c_v[1]
    for i in 1 : n
        bk[i] = ce + x_v[i] * bk1_v[i] - bk2_v[i]
    end
    
    bk
    
#     ## Hack to corect the fact bk_v keeps changing    
#     for i in 1 : n
#         bk_v[i] = ce + .5 * x_v[i] * bk1_v[i] - bk2_v[i]
#     end
#     
# 
#      j = mod(length(c), 3)
#     j == 0 ? bk1 : j == 1 ? bk : bk2
end