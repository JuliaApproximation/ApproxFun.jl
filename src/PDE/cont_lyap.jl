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


# function cont_reduce_dofs{T<:Fun}( R,G::Vector{T}, A::Array, M::Operator, F::Array )
#     if length(R) > 0
#         # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
#         # then kill the row by subtracting
#         # MXR'[:,k]*A'[k,:]  from MXA'
#         # i.e., subtacting A[:,k]*R[k,:] from A
#         # and M*G'[:,k]*A'[k,:] from F
#         # i.e. M*G[k]*A[:,k]' from 
#         
#         for k = 1:size(R,1)
#             MG = M*G[k].coefficients         # coefficients in the range space of M      
#             MGA = MG*A[:,k].'
#             m=max(size(F,1),size(MGA,1))
#             F = pad(F,m,size(F,2)) - pad(MGA,m,size(F,2))
#             A = A - A[:,k]*R[k,:]
#         end
#     end
#         
#     A, F
# end

function adaptiveminus(F::Matrix,G::Matrix)
    m=max(size(F,1),size(G,1))
    n=max(size(F,2),size(G,2))
    pad(F,m,n) - pad(G,m,n)    
end

function cont_reduce_dofs{T<:Fun,NT<:Number}( A::AbstractArray{NT},M::Operator,G::Vector{T},F::AbstractArray )
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from 
        
    for k = 1:length(G)
        MG = M*G[k].coefficients         # coefficients in the range space of M      
        MGA = MG*full(A[:,k]).'
        F = adaptiveminus(F,MGA)
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
cont_constrained_lyapuptriang{OSS,T,FT}(OS::PDEOperatorSchur{OSS,T},Gx,F::Array{FT},nx=100000)=cont_constrained_lyapuptriang(promote_type(T,FT),OS,Gx,F,nx)
#cont_constrained_lyapuptriang{N}(::Type{N},OS::PDEOperatorSchur,Gx,F::Array)=cont_constrained_lyapuptriang(N,OS,Gx,F,100000)


function cont_constrained_lyap{OSS<:DiagonalOperatorSchur,T}(OS::PDEOperatorSchur{OSS},Gxin,Gyin,F::Matrix{T},nx=100000)    
    n = size(OS.S,1)    
    F=pad(F,size(F,1),n)
    Gx=toarray(Gxin,n)    
    
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
    Gx=toarray(Gxin,n)    
    TYP=promote_type(eltype(OS),T)
    Y=Array(Fun{typeof(domainspace(OS.Rdiags[1])),TYP},n) 


    for k=1:n
        op=OS.Rdiags[k]
        rhs=T[Gx[:,k]...,F[:,k]...]
        Y[k]=chop!(linsolve([OS.Bx[k],op],rhs;maxlength=nx),eps())
    end  
    
    Y   
end



function cont_constrained_lyapuptriang{N,OSS<:OperatorSchur}(::Type{N},OS::PDEOperatorSchur{OSS},Gx,F::Array,nx::Integer)
    n = min(size(OS.S.T,2),max(size(F,2),size(Gx,2)))
    F=pad(F,size(F,1),n)
    if !isempty(Gx)
        Gx=pad(Gx,size(Gx,1),n)
    end


    Y=Array(Fun{typeof(domainspace(OS,1)),N},n)
    PY=Array(Vector{N},n)
    SY=Array(Vector{N},n)

    k=n
    m=n  # max length
    
    
    while k>1
        if OS.S.T[k,k-1] == 0 && OS.S.R[k,k-1] == 0        
            rhs = pad!(F[:,k],m)
            if k < n
                for j=k+1:n
                    for s=1:length(PY[j])
                        rhs[s]-=OS.S.R[k,j]*PY[j][s]
                    end
                    for s=1:length(SY[j])
                        rhs[s]-=OS.S.T[k,j]*SY[j][s]
                    end                    
                end
            end

            op=OS.Rdiags[k]
            if isempty(Gx)
                Y[k]=chop!(linsolve([OS.Bx,op],rhs;maxlength=nx),eps())            
            else
                Y[k]=chop!(linsolve([OS.Bx,op],[Gx[:,k],rhs];maxlength=nx),eps())
            end
            
            PY[k]=OS.Lx*Y[k].coefficients;SY[k]=OS.Mx*Y[k].coefficients
            m=max(m,length(PY[k]),length(SY[k]))
            
            k-=1
        else
            rhs1=F[:,k-1]
            rhs2=F[:,k]
        
            if k < n
                for j=k+1:n
                    rhs1= adaptiveminus!(rhs1,OS.S.R[k-1,j]*PY[j],OS.S.T[k-1,j]*SY[j])
                    rhs2= adaptiveminus!(rhs2,OS.S.R[k,j]*PY[j],OS.S.T[k,j]*SY[j])        
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
        
            PY[k-1]=OS.Lx*Y[k-1].coefficients; PY[k]=OS.Lx*Y[k].coefficients
            SY[k-1]=OS.Mx*Y[k-1].coefficients; SY[k]=OS.Mx*Y[k].coefficients
            
            m=max(m,length(PY[k]),length(SY[k]),length(PY[k-1]),length(SY[k-1]))  
            
            k-=2
        end
    end

    if k == 1
        rhs = F[:,k]

        for j=2:n
            rhs= adaptiveminus!(rhs,OS.S.R[k,j]*PY[j],OS.S.T[k,j]*SY[j])
        end

        Y[k]=chop!([OS.Bx,OS.Rdiags[k]]\[Gx[:,k],rhs],eps())
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




function cont_constrained_lyap{OSS<:OperatorSchur}(OS::PDEOperatorSchur{OSS},Gx,Gyin,F::Array,nx=100000)    
    Gy=regularize_bcs(OS.S,Gyin)
    F=cont_reduce_dofs(OS.S,OS.Lx,OS.Mx,Gy,F)  
    
    Q2 = OS.S.Q
    
    F=pad(F,size(F,1),size(Q2,1))*Q2
    
    
    ny=size(OS.S,2)
    Ky=numbcs(OS.S)
    
    ## we've discretized, in y, and rhs for Bx is a function of y
    # so we need to discetize it as well
    # and permute columns by P

    if !isempty(Gx)
        Gx=toarray(Gx,ny)*OS.S.bcP  
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



# function cont_constrained_lyap(Bxin,Byin,Lin,Min,F::Array,ny)
#     Xop=promotespaces([Lin[1],Min[1]])
#     Lx=SavedBandedOperator(Xop[1]);Mx=SavedBandedOperator(Xop[2])
# 
#     #discretize in Y
#     By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny) 
#     Ry,Gy,Ly,My,Py=regularize_bcs(By,Gy,full(Ly),full(My)) 
#     Ly,F = cont_reduce_dofs(Ry,Gy,Ly,Lx,F)     
#     My,f = cont_reduce_dofs(Ry,Gy,My,Mx,F)     
# 
#     Ky = size(By,1)
#     B=Ly[:,Ky+1:end]
#     D=My[:,Ky+1:end]
#     BD=schurfact(full(B),full(D))
#     Q2=BD[:left];Z2=BD[:right]
#     R=BD[:S]; T=BD[:T]
# 
#     F=pad(F,size(F,1),size(Q2,1))*Q2
#     
# 
#     ## we've discretized, in y, and rhs for Bx is a function of y
#     # so we need to discetize it as well
#     Gx=toarray(Bxin[2],ny)
#     # remove unused DOFs and rearrange columns
#     Gx=Gx[:,Ky+1:end]*Z2
# 
# 
#     Y=cont_constrained_lyapuptriang(Bxin[1],Gx,Lx,R,Mx,T,F)
#     
#     X22=Z2*Y  #think of it as transpose
#     X11=convert(typeof(X22),Gy-Ry[:,Ky+1:end]*X22) #temporary bugfix since Gy might have worse type
#     [X11,X22].'
# 
# end


