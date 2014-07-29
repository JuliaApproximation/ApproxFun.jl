export pdesolve


include("cont_lyap.jl")


##TODO: Unify with coefficients()
toarray{T<:Functional}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
toarray{T<:IFun}(B::Array{T},n)=Float64[    j<=length(B[k])?B[k].coefficients[j]:0 for  k=1:length(B),j=1:n];
function toarray(B::Array,n)
    ret = zeros(length(B),n)
    
    for k=1:length(B), j=1:n
        if  typeof(B[k]) <: IFun
            ret[k,j] = j<=length(B[k])?B[k].coefficients[j]:0 
        elseif typeof(B[k]) <: Number && j == 1
            ret[k,j] = B[k]
        end
    end

    ret
end

function toarray{T<:Operator}(A::Vector{T},n::Integer,m::Integer)
    ret = zeros(n,m)
    
    nbc = typeof(A[end])<:Functional?length(A):length(A)-1
    for k=1:nbc
        ret[k,:]=A[k][1:m]
    end
    
    if nbc < length(A)
        ret[nbc+1:end,:]=A[end][1:n-nbc,1:m]
    end
    
    ret
end

function pdetoarray(Byin,Lin,Min,ny::Integer)
    Yop=promotespaces([Lin,Min])

    
    By=toarray(Byin[1],ny)   
    nbcy=length(Byin[1])
    
    Ly=Yop[1][1:ny-nbcy,1:ny]
    My=Yop[2][1:ny-nbcy,1:ny]    
    
    By,Byin[2],Ly,My    
end

function pdetoarray(Bxin,Byin,Lin,Min,nx::Integer,ny::Integer)
    Bx,Gx,Lx,Mx=pdetoarray(Bxin,Lin[1],Min[1],nx)
    By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny)    
    

    Gx=toarray(Gx,ny);Gy=toarray(Gy,nx)    
    
    
    Bx,Gx,By,Gy,Lx,Ly,Mx,My    
end




function pdesolve(Bxin,Byin,Lin,Min,Fin::Number,nx::Integer,ny::Integer)
    F=zeros(nx-length(Bxin[1]),ny-length(Byin[1]))
    F[1,1]=Fin
    
    pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
end

function pdesolve(Bxin,Byin,Lin,Min,Fin::Fun2D,nx::Integer,ny::Integer)
    Xop=promotespaces([Lin[1],Min[1]])
    Yop=promotespaces([Lin[2],Min[2]])
    
    Xsp=rangespace(Xop[1])
    Ysp=rangespace(Yop[1])    
    
    nbcx=length(Bxin[1]);nbcy=length(Byin[1])    
    
    F=pad(coefficients(Fin,Xsp,Ysp),nx-nbcx,ny-nbcy)
    pdesolve(Bxin,Byin,Lin,Min,F,nx,ny)
end



function pdesolve(Bxin,Byin,Lin,Min,F::Array,nx::Integer,ny::Integer)

    Bx,Gx,By,Gy,Lx,Ly,Mx,My=pdetoarray(Bxin,Byin,Lin,Min,nx,ny)
    
    
    Fun2D(constrained_lyap({Bx Gx; By Gy},{Lx,Ly},{Mx,My},F),
            domain(Lin[1]),domain(Lin[2]))
end


function pdesolve(Bxin,Byin,Lin,Min,Fin::Number,ny::Integer)
    F=zeros(1,ny-length(Byin[1]))
    F[1,1]=Fin
    
    pdesolve(Bxin,Byin,Lin,Min,F,ny)
end

pdesolve(Bxin,Byin,Lin,Min,Fin::Fun2D,ny::Integer)=pdesolve(Bxin,Byin,Lin,Min,coefficients(Fin,Xsp,Ysp),ny)



pdesolve(Bx,By,Lin,Min,F::Array,ny::Integer)=Fun2D(cont_constrained_lyap(Bx, By,Lin,Min,F,ny),domain(Lin[2]))


