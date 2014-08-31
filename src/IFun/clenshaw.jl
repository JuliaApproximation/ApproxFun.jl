
type ClenshawPlan{T}
    bk::Vector{T}
    bk1::Vector{T}
    bk2::Vector{T}
end

ClenshawPlan{T}(::Type{T},n::Integer)=ClenshawPlan(Array(T,n),Array(T,n),Array(T,n))


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
    clenshaw(c,x,ClenshawPlan(T,n))
end


#Clenshaw routine for many Funs, x is a vector of same number of funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Vector{T})=clenshaw(c,x,ClenshawPlan(T,size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Vector{T},plan::ClenshawPlan{T})
    bk=plan.bk    
    bk1=plan.bk1
    bk2=plan.bk2


    n=size(c)[2] #number of funs
    m=size(c)[1] #number of coefficients

    
    for i = 1:n
        @inbounds bk1[i] = 0.
        @inbounds bk2[i] = 0.
        @inbounds bk[i] = 0.        
    end

    for k=m:-1:2
        for j=1:n
            ck = c[k,j]

            @inbounds bk[j] = ck + 2x[j] * bk1[j] - bk2[j]
        end
        
        bk2, bk1, bk = bk1, bk, bk2
    end


    for j = 1:n
        ce = c[1,j]
        @inbounds bk[j] = ce + x[j] * bk1[j] - bk2[j]
    end
    
    bk    
end



#Clenshaw routine for many Funs
#each fun is a column
clenshaw{T<:Number}(c::Array{T,2},x::Real)=clenshaw(c,x,ClenshawPlan(T,size(c)[2]))
function clenshaw{T<:Number}(c::Array{T,2},x::Real,plan::ClenshawPlan{T})
    bk=plan.bk    
    bk1=plan.bk1
    bk2=plan.bk2


    n=size(c)[2] #number of funs
    m=size(c)[1] #number of coefficients

    for i = 1:n
        @inbounds bk1[i] = 0.
        @inbounds bk2[i] = 0.
        @inbounds bk[i] = 0.                
    end

    for k=m:-1:2
        for j=1:n
            @inbounds ck = c[k,j]

            @inbounds bk[j] = ck + 2x * bk1[j] - bk2[j]
        end
        
        bk2, bk1, bk = bk1, bk, bk2
    end


    for j = 1:n
        ce = c[1,j]
        @inbounds bk[j] = ce + x * bk1[j] - bk2[j]
    end
    
    bk    
end


#Clenshaw routine for many points
#Note that bk1, bk2, and bk are overwritten
function clenshaw{T<:Number,M<:Real}(c::Vector{T},x::Vector{M},plan::ClenshawPlan{T})
    bk=plan.bk    
    bk1=plan.bk1
    bk2=plan.bk2

    n = length(x)
#    x=2x

    
    for i = 1:n
        @inbounds bk1[i] = 0.
        @inbounds bk2[i] = 0.
        @inbounds bk[i] = 0.                
    end

    for k in  length(c):-1:2
        ck = c[k]
        for i in 1 : n
            @inbounds bk[i] = ck + 2x[i] * bk1[i] - bk2[i]
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ce = c[1]
    for i in 1 : n
        @inbounds  bk[i] = ce + x[i] * bk1[i] - bk2[i]
    end
    
    bk
end


# overwrite x
clenshaw!{T<:Real}(c::Vector{T},x::Vector{T})=clenshaw!(c,x,ClenshawPlan(T,length(x)))
function clenshaw!{T<:Real}(c::Vector{T},x::Vector{T},plan::ClenshawPlan{T})
    bk=plan.bk    
    bk1=plan.bk1
    bk2=plan.bk2

    n = length(x)
#    x=2x

    
    for i = 1:n
        @inbounds bk1[i] = 0.
        @inbounds bk2[i] = 0.
        @inbounds bk[i] = 0.                
    end

    for k in  length(c):-1:2
        ck = c[k]
        for i in 1 : n
            @inbounds bk[i] = ck + 2x[i] * bk1[i] - bk2[i]
        end
        bk2, bk1, bk = bk1, bk, bk2
    end

    ce = c[1]
    for i in 1 : n
        @inbounds  x[i] = ce + x[i] * bk1[i] - bk2[i]
    end
    
    x
end

