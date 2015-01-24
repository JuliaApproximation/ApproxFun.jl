function adaptiveplus(f::Vector,g::Vector)
    if length(f)>length(g)
        ret=copy(f)
        @inbounds ret[1:length(g)]+=g
    else
        ret=copy(g)
        @inbounds ret[1:length(f)]+=f
    end
    
    ret
end

##May modify either f or g
function adaptiveplus!(f::Vector,g::Vector)
    if length(f)>length(g)
        @simd for k=1:length(g)
            @inbounds f[k]+=g[k]
        end
        f
    else
        @simd for k=1:length(f)
            @inbounds g[k]+=f[k]
        end
        g
    end
end

function adaptiveminus!(f::Vector,g::Vector)
    if length(f)>length(g)
        @simd for k=1:length(g)
            @inbounds f[k]-=g[k]
        end
        f
    else
        @simd for k=1:length(g)
            g[k]*=-1
        end
        @simd for k=1:length(f)
            @inbounds g[k]+=f[k]
        end
        g
    end
end
        
##May modify either f,g or h
adaptiveplus!(f,g,h)=adaptiveplus!(adaptiveplus!(f,g),h)

#f-g-h
adaptiveminus!(f,g,h)=adaptiveminus!(adaptiveminus!(f,g),h)


#Use XR' = G' = [G1 G2 G3...] to reduce columns of A in
# MXA' + *X* =F
# here G is a vector of Funs


function adaptiveminus(F::Matrix,G::Matrix)
    m=max(size(F,1),size(G,1))
    n=max(size(F,2),size(G,2))
    pad(F,m,n) - pad(G,m,n)    
end

function cont_reduce_dofs!{T<:Fun,NT<:Number}( A::AbstractArray{NT},M::Operator,G::Vector{T},F::ProductFun )
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from 
        
    for k = 1:length(G)
        MG = M*G[k]         # coefficients in the range space of M      
        for j=1:length(F.coefficients)
            F.coefficients[j]-=A[j,k]*MG
        end
    end
        
    F
end


# G is  ∞ x K array
# A is ∞ x K list of opcols
# M is ∞ x ∞ operator
function cont_reduce_dofs{NT<:Number,T<:Number}( A::AbstractArray{NT},M::AbstractArray{T},G::Array,F::AbstractArray )
    MGA=M*pad(G,size(M,1),size(G,2))*full(A).'
    adaptiveminus(F,MGA)
end


function cont_reduce_dofs{T<:Fun,NT<:Number,MT<:Number}( A::AbstractArray{NT},M::AbstractArray{MT},G::Vector{T},F::AbstractArray )
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from 
        
    for k = 1:length(G)
        MG = M*pad(G[k].coefficients,size(M,2))         # coefficients in the range space of M      
        MGA = MG*full(A[:,k]).'
        F=adaptiveminus(F,MGA)
    end
        
    F
end

function cont_reduce_dofs{T<:Fun}(S::OperatorSchur,L::Operator,M::Operator,G::Vector{T},F::AbstractArray)
    F=cont_reduce_dofs(S.Lcols,L,G,F)
    cont_reduce_dofs(S.Mcols,M,G,F)    
end

function cont_reduce_dofs!{T<:Fun}(S::OperatorSchur,L::Operator,M::Operator,G::Vector{T},F::ProductFun)
    F=cont_reduce_dofs!(S.Lcols,L,G,F)
    cont_reduce_dofs!(S.Mcols,M,G,F)    
end



function cont_reduce_dofs{M<:AbstractArray}(Ax::Vector{M},Ay::Vector,G,F::AbstractArray)
    @assert length(Ax)==length(Ay)
    for k=1:length(Ax)
        F=cont_reduce_dofs(Ax[k],Ay[k],G,F)
    end
    
    F
end




regularize_bcs(S::OperatorSchur,Gy)=length(Gy)==0?Gy:S.bcQ*Gy


# Solve Bx*Y=Gx and P*Y*R' + S*Y*T' = F 
# where R and T are upper triangular

##TODO: Support complex in boundary conditions
cont_constrained_lyapuptriang{OSS,T}(OS::PDEOperatorSchur{OSS,T},Gx,F::ProductFun,nx=100000)=cont_constrained_lyapuptriang(promote_type(T,eltype(F)),OS,Gx,F,nx)
#cont_constrained_lyapuptriang{N}(::Type{N},OS::PDEOperatorSchur,Gx,F::Array)=cont_constrained_lyapuptriang(N,OS,Gx,F,100000)


function cont_constrained_lyap{OSS<:DiagonalOperatorSchur,T}(OS::PDEOperatorSchur{OSS},Gxin,Gyin,F::Matrix{T},nx=100000)    
    n = size(OS.S,1)    
    F=pad(F,size(F,1),n)
    Gx=pad(coefficients(Gxin).',:,n)
    
    TYP=promote_type(eltype(OS),T)    
    Y=Array(Fun{typeof(domainspace(OS,1)),TYP},n)


    for k=1:n
        op=OS.Rdiags[k]
        rhs=T[Gx[:,k]...,F[:,k]...]
        Y[k]=chop!(linsolve([OS.Bx,op],rhs;maxlength=nx),eps())
    end  
    
    Y   
end

function cont_constrained_lyap{T}(OS::PDEProductOperatorSchur,Gxin,Gyin,F::Matrix{T},nx=100000)    
    n = length(OS.Rdiags)
    F=pad(F,size(F,1),n)
    Gx=pad(coefficients(Gxin).',:,n)
    TYP=promote_type(eltype(OS),T)
    Y=Array(Fun{typeof(domainspace(OS.Rdiags[1])),TYP},n) 


    for k=1:n
        op=OS.Rdiags[k]
        rhs=T[Gx[:,k]...,F[:,k]...]
        Y[k]=chop!(linsolve([OS.Bx[k],op],rhs;maxlength=nx),eps())
    end  
    
    Y   
end



function cont_constrained_lyapuptriang{N,OSS<:OperatorSchur}(::Type{N},OS::PDEOperatorSchur{OSS},Gx,F::ProductFun,nx::Integer)
    n = min(size(OS.S.T,2),max(size(F,2),size(Gx,2)))

    if !isempty(Gx)
        Gx=pad(Gx,size(Gx,1),n)
    end

    rs=rangespace(OS)
    @assert space(F)==rs

    Y=Array(Fun{typeof(domainspace(OS,1)),N},n) # the solution
    PY=Array(Fun{typeof(rs[1]),N},n) # the x-op times the solution
    SY=Array(Fun{typeof(rs[2]),N},n) # the y-op times the solution

    k=n
    m=n  # max length
    
    
    while k≥1
        if k==1 || (OS.S.T[k,k-1] == 0 && OS.S.R[k,k-1] == 0        )  # triangular setting
            rhs = k≤length(F.coefficients)?F.coefficients[k]:zeros(rs[1])
            
            if k < n
                for j=k+1:n
                    axpy!(-OS.S.R[k,j],PY[j],rhs)
                    axpy!(-OS.S.T[k,j],SY[j],rhs)
                end
            end

            op=OS.Rdiags[k]
            if isempty(Gx)
                Y[k]=chop!(linsolve([OS.Bx,op],rhs;maxlength=nx),eps())            
            else
                Y[k]=chop!(linsolve([OS.Bx,op],[Gx[:,k],rhs];maxlength=nx),eps())
            end
            
            if k > 1
                PY[k]=OS.Lx*Y[k];SY[k]=OS.Mx*Y[k]
            end
            
            k-=1
        else # quasitriangular
            rhs1 = k-1≤length(F.coefficients)?F.coefficients[k-1]:zeros(rs[1])        
            rhs2 = k≤length(F.coefficients)?F.coefficients[k]:zeros(rs[1])                    
        
            if k < n
                for j=k+1:n
                    axpy!(-OS.S.R[k-1,j],PY[j],rhs1)
                    axpy!(-OS.S.T[k-1,j],SY[j],rhs1)         

                    axpy!(-OS.S.R[k,j],PY[j],rhs2)
                    axpy!(-OS.S.T[k,j],SY[j],rhs2)                                        
                end
            end
        
        
            A=[blkdiag(OS.Bx,OS.Bx);
                OS.S.R[k-1:k,k-1:k].*OS.Lx+OS.S.T[k-1:k,k-1:k].*OS.Mx]
            if isempty(Gx)
                b=Any[rhs1,rhs2]
            else
                b=Any[Gx[:,k-1]...,Gx[:,k]...,rhs1,rhs2]
            end
            y=vec(linsolve(A,b;maxlength=nx))
            Y[k-1]=chop!(y[1],eps());Y[k]=chop!(y[2],eps())
        
            PY[k-1]=OS.Lx*Y[k-1]; PY[k]=OS.Lx*Y[k]
            SY[k-1]=OS.Mx*Y[k-1]; SY[k]=OS.Mx*Y[k]
            
            k-=2
        end
    end
    
    Y
end


## Solve ∞-dimensional lyap equation
#
#   Bx*X=[Gx[1](y);Gx[2](y);...]=Gx'
#   X*By' = [Gy[1](x) Gy[2](x) ...]
#   Lx*X*Ly' + Mx*X*My' = F
#
# by discretizing in y

#TODO: F should be adaptive rather than array
#TODO: Describe precisely the permutation structure of OS so this algorithm actually makes sense
#      to others (and me too!!)




function cont_constrained_lyap{OSS<:OperatorSchur}(OS::PDEOperatorSchur{OSS},Gx::Vector,Gyin::Vector,F::ProductFun,nx=100000)    
    Gy=regularize_bcs(OS.S,Gyin)
    F=cont_reduce_dofs!(OS.S,OS.Lx,OS.Mx,Gy,F)  
    
     # Q2 says how to rearrange the columns of F so that the operator is upper triangular
    Q2 = OS.S.Q 
    F=ProductFun(Q2[1:length(F.coefficients),:].'*F.coefficients,space(F))
    
    
    ny=size(OS.S,2)
    Ky=numbcs(OS.S)
    
    ## we've discretized, in y, and rhs for Bx is a function of y
    # so we need to discetize it as well
    # and permute columns by P

    if !isempty(Gx)
        # bcP recombines boundary conditions        
        Gx=pad(coefficients(Gx).',:,ny)*OS.S.bcP  
        # remove unused DOFs and rearrange columns
        Gx=Gx[:,Ky+1:end]*OS.S.Z
    else
        Gx=[]
    end
    
    Y=cont_constrained_lyapuptriang(OS,Gx,F,nx)
    # Y is a Vector{Fun}, so that Y[k][j] corresponds to matrix element M[k,j]
    # This means acting on Y is acting on *columns* of M
    
    X22=OS.S.Z*Y  #think of it as transpose
    
    X11=convert(typeof(X22),Gy-OS.S.bcs[:,Ky+1:end]*X22) #temporary bugfix since Gy might have worse type
    X=[X11,X22]    
    X=OS.S.bcP*X        # this is equivalent to acting on columns by P'
end

