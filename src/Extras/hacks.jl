## Constructors that involve MultivariateFun

function Fun(f::Function,d::BivariateSpace)
    if f==zero
        zeros(d)
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        Fun(ProductFun(f,d))
    else
        Fun(ProductFun((x,y)->f((x,y)),d))
    end
end


function Fun(f::Function,d::BivariateSpace,n::Integer)
    if (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        defaultFun(x->f(x...),d,n)
    else
        defaultFun(f,d,n)
    end
end

function Fun(f::Function)
    if (isgeneric(f)&&applicable(f,0)) || (!isgeneric(f)&&arglength(f)==1)
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
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
            Fun(f,Interval()^2)
    else
        error("Function not defined on interval or square")
    end
end






## ConstatnOPerators can always be promoted
promotedomainspace{T,CO<:ConstantOperator}(S::SpaceOperator{T,CO},sp::AnySpace)=S
promotedomainspace{T,CO<:ConstantOperator}(S::SpaceOperator{T,CO},sp::UnsetSpace)=S
promotedomainspace{T,CO<:ConstantOperator}(S::SpaceOperator{T,CO},sp::Space)=SpaceOperator(S.op,sp,sp)





## MutableOperator bor BandedMatrix

function resizedata!{T<:Matrix,M<:BandedOperator,R}(B::MutableOperator{T,M,R},n::Integer)
    resizedata!(B.fill,n)

    if n > B.datalength
        nbc=B.fill.numbcs

        if n > size(B.data,1)
            newdata=blockbandzeros(eltype(T),n,:,bandinds(B.data),blockbandinds(B.op))

            for k=1:B.datalength,j=max(1,k+bandinds(B.data,1)):k+bandinds(B.data,2)
                newdata.data[k,j]=B.data.data[k,j]
            end
            B.data=newdata
        end

        addentries!(B.op,IndexStride(B.data,nbc,0),B.datalength+1-nbc:n-nbc,:)
        B.datalength = n
    end

    B
end

function givensmatrix(a::Matrix,b::Matrix)
    q,r=qr([a;b];thin=false)
    q[1:size(a,1),1:size(a,2)],q[1:size(a,1),size(a,1)+1:end],q[size(a,1)+1:end,1:size(a,2)],q[size(a,1)+1:end,size(a,2)+1:end]
end


#TOD: Bcs
function unsafe_getindex{T<:Matrix,R}(B::FillMatrix{T,R},k::Integer,j::Integer)
    zeros(eltype(T),k,j)
end


#TODO: Fix hack override
function pad{T<:Vector}(f::Vector{T},n::Integer)
	if n > length(f)
	   ret=Array(T,n)
	   ret[1:length(f)]=f
	   for j=length(f)+1:n
	       ret[j]=zeros(eltype(T),j)
	   end
       ret
	else
        f[1:n]
	end
end




function backsubstitution!{T<:Vector}(B::MutableOperator,u::Array{T})
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.fill.numbcs
    A=B.data

    @assert nbc==0

    for c=1:size(u,2)

        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            @simd for j=k+1:n
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.l+1,k]\u[k,c]
        end

       #filled rows
        for k=n-b-1:-1:1
            @simd for j=k+1:k+b
                @inbounds u[k,c]-=A.data[j-k+A.l+1,k]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.l+1,k]\u[k,c]
        end
    end
    u
end



function applygivens!(ca::Matrix,cb,mb,a,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    @simd for j = jr
        @inbounds B1 = B.data[j-k1+B.l+1,k1]    #B[k1,j]
        @inbounds B2 = B.data[j-k2+B.l+1,k2]    #B[k2,j]

        @inbounds B.data[j-k1+B.l+1,k1]=ca*B1
        @inbounds B.data[j-k2+B.l+1,k2]=a*B2
        BLAS.gemm!('N','N',1.0,cb,B2,1.0,B.data[j-k1+B.l+1,k1])
        BLAS.gemm!('N','N',1.0,mb,B1,1.0,B.data[j-k2+B.l+1,k2])
    end

    B
end

function applygivens!(ca::Matrix,cb,mb,a,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        @inbounds B2 = B.data[j-k2+B.l+1,k2]   #B[k2,j]
        @inbounds B.data[j-k2+B.l+1,k2]=a*B2
    end

    B
end

function applygivens!(ca::Matrix,cb,mb,a,B::Matrix,k1::Integer,k2::Integer)
    for j = 1:size(B,2)
        @inbounds B1 = B[k1,j]
        @inbounds B2 = B[k2,j]

        @inbounds B[k1,j],B[k2,j]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end


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
