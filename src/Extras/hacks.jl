## Functions that depend on the structure of BandedMatrix


function pad!(A::BandedMatrix,n,m)
    A.data=pad(A.data,size(A.data,1),m)
    A.m=n
    A
end



# default copy is to loop through
# override this for most operators.
function Base.copy(S::SubBandedOperator)
    l,u = bandwidth(S,1),bandwidth(S,2)
    Y=BandedMatrix(eltype(S),size(S,1),size(S,2),l,u)

    for (k,j) in eachbandedindex(S)
        @inbounds Y.data[k-j+u+1,j]=X[k,j]
    end

    Y
end


# linear algebra


function toeplitz_addentries!(v::Vector,A::BandedMatrix,kr::UnitRange)
    if !isempty(v)
        @inbounds v1=v[1]
        @simd for k=kr
            @inbounds A.data[A.u+1,k]+=2v1
        end

        for j=2:length(v)
            @inbounds vj=v[j]
            @simd for k = kr
                @inbounds A.data[A.u-j+2,k+j-1] +=vj
            end
            @simd for k = max(kr[1],j):kr[end]
                @inbounds A.data[A.u+j,k-j+1]+=vj
            end
        end
    end
    A
end

# ca gives conj(a), cb gives conj(b), mb gives -b
# This is written to support matrix value and a and b

# applygivens! considers three cases:
#  The case where two entries in the banded part are manipulated
#  The case where one enty in the band part and one in the filled in part are manipulated
#  The case where the fill is manipulated
#

function applygivens!(ca,cb,mb,a,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    @simd for j = jr
        @inbounds B1 = B.data[k1-j+B.u+1,j]    #B[k1,j]
        @inbounds B2 = B.data[k2-j+B.u+1,j]    #B[k2,j]

        @inbounds B.data[k1-j+B.u+1,j],B.data[k2-j+B.u+1,j]= ca*B1 + cb*B2,mb*B1 + a*B2
    end

    B
end

function applygivens!(ca,cb,mb,a,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        B1 = unsafe_getindex(F,k1,j)
        @inbounds B2 = B.data[k2-j+B.u+1,j]   #B[k2,j]

        @inbounds B.data[k2-j+B.u+1,j]=mb*B1 + a*B2
    end

    B
end

function givensreduceab!{T,M,R}(B::MutableOperator{T,M,R},k1::Integer,k2::Integer,j1::Integer)
    bnd=B.bandinds
    A=B.data

    @inbounds a=A.data[k1-j1+A.u+1,j1]  #A[k1,j1]
    @inbounds b=A.data[k2-j1+A.u+1,j1]  #A[k2,j1]

    ca,cb,mb,a=givensmatrix(a,b)


    #Assuming that left rows are already zero
    applygivens!(ca,cb,mb,a,B.data,k1,k2,j1:k1+B.bandinds[end])
    applygivens!(ca,cb,mb,a,B.fill,B.data,k1,k2,k1+B.bandinds[end]+1:k2+B.bandinds[end])
    applygivens!(ca,cb,mb,a,B.fill.data,k1,k2)


    ca,cb,mb,a
end


function backsubstitution!(B::MutableOperator,u::Array)
    n=size(u,1)
    b=B.bandinds[end]
    nbc = B.fill.numbcs
    A=B.data
    T=eltype(u)
    pk = zeros(T,nbc)

    for c=1:size(u,2)
        fill!(pk,zero(T))

        # before we get to filled rows
        for k=n:-1:max(1,n-b)
            @simd for j=k+1:n
                @inbounds u[k,c]-=A.data[k-j+A.u+1,j]*u[j,c]
            end

            @inbounds u[k,c] /= A.data[A.u+1,k]
        end

       #filled rows
        for k=n-b-1:-1:1
            @simd for j=1:nbc
                @inbounds pk[j] += u[k+b+1,c]*B.fill.bc[j].data[k+b+1]
            end

            @simd for j=k+1:k+b
                @inbounds u[k,c]-=A.data[k-j+A.u+1,j]*u[j,c]
            end

            @simd for j=1:nbc
                @inbounds u[k,c] -= B.fill.data[k,j]*pk[j]
            end

            @inbounds u[k,c] /= A.data[A.u+1,k]
        end
    end
    u
end








function kronaddentries!(A,B,M,kr::Range)
    m=max(size(A,2),size(B,2))
    l=A.l+B.l;u=A.u+B.u

    for k=kr,j=max(1,k-l):k+u
        nl=min(A.l,B.u+k-j);nu=min(A.u,B.l+j-k)
        @inbounds Mkj=M[k,j]
        for κ=1:k,ξ=max(1,κ-nl):min(j,κ+nu)
            #Mkj[κ,ξ]+=A[κ,ξ]*B[k-κ+1,j-ξ+1]
            @inbounds Mkj[κ,ξ]+=A.data[κ-ξ+A.u+1,ξ]*B.data[k-κ-j+ξ+B.u+1,j-ξ+1]
        end
    end
    M
end


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
    elseif (isgeneric(f)&&applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
            Fun(f,Interval()^2)
    else
        error("Function not defined on interval or square")
    end
end







## ConstatnOPerators can always be promoted
promotedomainspace{CO<:ConstantOperator}(S::SpaceOperator{CO},sp::AnySpace)=S
promotedomainspace{CO<:ConstantOperator}(S::SpaceOperator{CO},sp::UnsetSpace)=S
promotedomainspace{CO<:ConstantOperator}(S::SpaceOperator{CO},sp::Space)=SpaceOperator(S.op,sp,sp)





## These hacks support PDEs with block matrices


# this is a hack since it assumes the growth in the blocks
function bzeros{T}(::Type{Matrix{T}},n::Integer,m::Integer,l::Integer,u::Integer)
    ret=BandedMatrix(Matrix{T},n,m,l,u)
    for k=1:n,j=max(1,k-l):min(m,k+u)
        ret[k,j]=zeros(T,k,j)  #The ::Number works around an 0.4 bug
    end
    ret
end




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
                @inbounds u[k,c]-=A.data[k-j+A.u+1,j]*u[j,c]   #A[k,j]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.u+1,k]\u[k,c]          #A[k,k]/u[k,c]
        end

       #filled rows
        for k=n-b-1:-1:1
            @simd for j=k+1:k+b
                @inbounds u[k,c]-=A.data[k-j+A.u+1,j]*u[j,c]
            end

            @inbounds u[k,c] = A.data[A.u+1,k]\u[k,c]
        end
    end
    u
end



function applygivens!(ca::Matrix,cb,mb,a,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    @simd for j = jr
        @inbounds B1 = B.data[k1-j+B.u+1,j]    #B[k1,j]
        @inbounds B2 = B.data[k2-j+B.u+1,j]    #B[k2,j]

        @inbounds B.data[k1-j+B.u+1,j]=ca*B1
        @inbounds B.data[k2-j+B.u+1,j]=a*B2
        BLAS.gemm!('N','N',1.0,cb,B2,1.0,B.data[k1-j+B.u+1,j])
        BLAS.gemm!('N','N',1.0,mb,B1,1.0,B.data[k2-j+B.u+1,j])
    end

    B
end

function applygivens!(ca::Matrix,cb,mb,a,F::FillMatrix,B::BandedMatrix,k1::Integer,k2::Integer,jr::Range)
    for j = jr
        @inbounds B2 = B.data[k2-j+B.u+1,j]   #B[k2,j]
        @inbounds B.data[k2-j+B.u+1,j]=a*B2
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
