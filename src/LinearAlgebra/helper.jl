## Helper routines
alternatingvector(n::Integer) = 2*mod([1:n],2) .- 1

function alternatesign!(v::Vector)
    n=length(v)
    for k=2:2:n
        v[k]=-v[k]
    end
    
    v
end

alternatesign(v::Vector)=alternatesign!(copy(v))



function pad!{T}(f::Vector{T},n::Integer)
	if (n > length(f))
		append!(f,zeros(T,n - length(f)))
	else
		resize!(f,n)
	end
end


function pad{T}(f::Vector{T},n::Integer)
	if (n > length(f))
        [f,zeros(T,n - length(f))]
	else
        f[1:n]
	end
end

function pad{T}(A::Array{T},n::Integer,m::Integer)
	if n <= size(A,1) && m <= size(A,2)
        A[1:n,1:m]
	else
        ret = zeros(T,n,m)
        
        if n <= size(A,1)
            ret[1:n,1:size(A,2)] = A[1:n,:]
        elseif m <= size(A,2)
            ret[1:size(A,1),1:m] = A[:,1:m]
        else
            ret[1:size(A,1),1:size(A,2)]=A
        end
        
        ret
	end
end




#TODO:padleft!

function padleft(f::Vector,n::Integer)
	if (n > length(f))
        [zeros(n - length(f)),f]
	else
        f[end-n+1:end]
	end
end



##chop!
function chop!(c::Vector,tol::Real)
    @assert tol > 0

    for k=[length(c):-1:1]
        if abs(c[k]) > tol
            resize!(c,k)
            return c
        end
    end
    
    resize!(c,0)
    c
end
