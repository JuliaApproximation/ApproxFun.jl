


#Use XR' = G' = [G1 G2 G3...] to reduce columns of A in
# MXA' + *X* =F
# here G is a vector of Funs

cont_reduce_dofs!{T<:Fun,NT<:Number}( A::AbstractArray{NT},M::Operator,G::Vector{T},F::Fun )=cont_reduce_dofs!( A,M,G,ProductFun(F))
function cont_reduce_dofs!{T<:Fun,NT<:Number}( A::AbstractArray{NT},M::Operator,G::Vector{T},F::ProductFun )
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from

    pad!(F,:,max(size(F,2),size(A,1)))

    for k = 1:length(G)
        MG = M*G[k]         # coefficients in the range space of M
        for j=1:min(length(F.coefficients),size(A,1))
            axpy!(-A[j,k],MG,F.coefficients[j]) # equivalent to X+=a*Y
        end
    end

    F
end


# G is  ∞ x K array
# A is ∞ x K list of opcols
# M is ∞ x ∞ operator
# used by kron
cont_reduce_dofs!{NT<:Number,T<:Number}( A::AbstractArray{NT},M::AbstractArray{T},G::Array,F::Fun )=cont_reduce_dofs!( A,M,G,ProductFun(F) )
function cont_reduce_dofs!{NT<:Number,T<:Number}( A::AbstractArray{NT},M::AbstractArray{T},G::Array,F::ProductFun )
    MGA=M*pad(G,size(M,1),size(G,2))*full(A).'
    pad!(F,:,max(size(F,2),size(MGA,2)))
    for j=1:size(MGA,2)
        axpy!(-1,MGA[:,j],F.coefficients[j])
    end
    F
end

cont_reduce_dofs!{T<:Fun,NT<:Number,MT<:Number}( A::AbstractArray{NT},M::AbstractArray{MT},G::Vector{T},F::Fun )=cont_reduce_dofs!(A,M,G,ProductFun(F))
function cont_reduce_dofs!{T<:Fun,NT<:Number,MT<:Number}( A::AbstractArray{NT},M::AbstractArray{MT},G::Vector{T},F::ProductFun )
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from

    for k = 1:length(G)
        MG = M*pad(G[k].coefficients,size(M,2))         # coefficients in the range space of M
        for j=1:length(F.coefficients)
            axpy!(-A[j,k],MG,F.coefficients[j]) # equivalent to X+=a*Y
        end
    end

    F
end


function cont_reduce_dofs!{T<:Fun}(S::OperatorSchur,L::Operator,M::Operator,G::Vector{T},F::ProductFun)
    F=cont_reduce_dofs!(S.Lcols,L,G,F)
    cont_reduce_dofs!(S.Mcols,M,G,F)
end

function cont_reduce_dofs!{T<:Fun}(S::OperatorSchur,L::Operator,M::Operator,G::Vector{T},F::Array)
    for Fk in F
        cont_reduce_dofs!(S,L,M,G,Fk)
    end
    F
end

function cont_reduce_dofs!{T<:Fun}(S::OperatorSchur,L::Operator,M::Operator,G::Matrix{T},F::Matrix)
    @assert size(G,2)==size(F,2)
    @assert size(F,1)==1


    for k=1:size(G,2)
        cont_reduce_dofs!(S,L,M,G[:,k],F[1,k])
    end
    F
end


cont_reduce_dofs!{M<:AbstractArray}(Ax::Vector{M},Ay::Vector,G,F::Fun)=cont_reduce_dofs!(Ax,Ay,G,ProductFun(F))
function cont_reduce_dofs!{M<:AbstractArray}(Ax::Vector{M},Ay::Vector,G,F::ProductFun)
    @assert length(Ax)==length(Ay)
    for k=1:length(Ax)
        F=cont_reduce_dofs!(Ax[k],Ay[k],G,F)
    end

    F
end





regularize_bcs(S::OperatorSchur,Gy) = length(Gy)==0?Gy:S.bcQ*Gy


# Solve Bx*Y=Gx and P*Y*R' + S*Y*T' = F
# where R and T are upper triangular

##TODO: Support complex in boundary conditions
cont_constrained_lyapuptriang{OSS,T}(OS::PDEOperatorSchur{OSS,T},Gx,F::ProductFun;kwds...)=cont_constrained_lyapuptriang(promote_type(T,eltype(F)),OS,Gx,F;kwds...)
#cont_constrained_lyapuptriang{N}(::Type{N},OS::PDEOperatorSchur,Gx,F::Array)=cont_constrained_lyapuptriang(N,OS,Gx,F,100000)

cont_constrained_lyap{OSS<:DiagonalOperatorSchur}(OS::PDEOperatorSchur{OSS},Gxin,Gyin,F::Fun;kwds...)=cont_constrained_lyap(OS,Gxin,Gyin,ProductFun(F);kwds...)
function cont_constrained_lyap{OSS<:DiagonalOperatorSchur}(OS::PDEOperatorSchur{OSS},Gxin,Gyin,F::ProductFun;kwds...)
    if isempty(OS.Bx)
        return cont_constrained_lyap_nobcs(OS,F;kwds...)
    end
    n = size(OS.S,1)
    F=pad(F,size(F,1),n)
    Gx=pad(coefficients(Gxin).',:,n)

    TYP=promote_type(eltype(OS),eltype(F))
    Y=Array(Fun{typeof(domainspace(OS,1)),TYP},n)


    for k=1:n
        op=OS.Rdiags[k]
        rhs=Any[Gx[:,k]...;F.coefficients[k]]
        Y[k]=chop!(linsolve(op,rhs;kwds...),eps())
    end

    Y
end

function cont_constrained_lyap_nobcs{OSS<:DiagonalOperatorSchur}(OS::PDEOperatorSchur{OSS},F::ProductFun;kwds...)
    n = size(OS.S,1)
    F=pad(F,size(F,1),n)

    TYP=promote_type(eltype(OS),eltype(F))
    Y=Array(Fun{typeof(domainspace(OS,1)),TYP},n)


    for k=1:n
        op=OS.Rdiags[k]
        rhs=F.coefficients[k]
        Y[k]=chop!(linsolve(op,rhs;kwds...),eps())
    end

    Y
end


function cont_constrained_lyapuptriang{N,OSS<:OperatorSchur}(::Type{N},OS::PDEOperatorSchur{OSS},Gx,F::ProductFun;kwds...)
    n = min(size(OS.S.T,2),max(size(F,2),size(Gx,2)))

    if !isempty(Gx)
        Gx=pad(Gx,size(Gx,1),n)
    end

    rs=rangespace(OS)
    @assert space(F)==rs

    Y=Array(Fun{typeof(domainspace(OS,1)),N},n) # the solution
    PY=Array(Fun{typeof(rs[1]),N},n) # the first x-op times the solution
    SY=Array(Fun{typeof(rs[1]),N},n) # the second x-op times the solution

    k=n
    m=n  # max length

    rhs=Array(Any,size(Gx,1)+1)
    TT=isempty(OS.Bx)?eltype(OS):promote_type(eltype(OS),mapreduce(eltype,promote_type,OS.Bx))

    blkops=Array(Operator{TT},2length(OS.Bx)+2,2)
    if !isempty(OS.Bx)
        blkops[1:2length(OS.Bx),:]=blkdiag(OS.Bx,OS.Bx)
    end
    blkrhs=Array(Any,2size(Gx,1)+2)


    while k≥1
        if k==1 || (OS.S.T[k,k-1] == 0 && OS.S.R[k,k-1] == 0        )  # triangular setting
            rhs[end] = k≤length(F.coefficients)?F.coefficients[k]:zeros(rs[1])

            if k < n
                for j=k+1:n
                    axpy!(-OS.S.R[k,j],PY[j],rhs[end]) # equivalent to X+=a*Y
                    axpy!(-OS.S.T[k,j],SY[j],rhs[end])
                end
            end



            for j=1:size(Gx,1)
                rhs[j]=Gx[j,k]
            end

            ops=OS.Rdiags[k]
            chop!(rhs[end],eps())
            Y[k]=chop!(linsolve(ops,rhs;kwds...),eps())

            if k > 1
                PY[k]=OS.Lx*Y[k];SY[k]=OS.Mx*Y[k]
            end

            k-=1
        else # quasitriangular
            rhs1 = k-1≤length(F.coefficients)?F.coefficients[k-1]:zeros(rs[1])
            rhs2 = k≤length(F.coefficients)?F.coefficients[k]:zeros(rs[1])

            if k < n
                for j=k+1:n
                    axpy!(-OS.S.R[k-1,j],PY[j],rhs1)  # +=
                    axpy!(-OS.S.T[k-1,j],SY[j],rhs1)  # +=

                    axpy!(-OS.S.R[k,j],PY[j],rhs2)  # +=
                    axpy!(-OS.S.T[k,j],SY[j],rhs2)  # +=
                end
            end

            blkops[end-1:end,:]=OS.S.R[k-1:k,k-1:k].*OS.Lx+OS.S.T[k-1:k,k-1:k].*OS.Mx

            for j=1:size(Gx,1)
                blkrhs[j]=Gx[j,k-1]
                blkrhs[j+size(Gx,1)]=Gx[j,k]
            end
            blkrhs[end-1]=rhs1;blkrhs[end]=rhs2


            y=vec(linsolve(blkops,blkrhs;kwds...))
            Y[k-1]=chop!(y[1],eps());Y[k]=chop!(y[2],eps())

            PY[k-1]=OS.Lx*Y[k-1]; PY[k]=OS.Lx*Y[k]
            SY[k-1]=OS.Mx*Y[k-1]; SY[k]=OS.Mx*Y[k]

            k-=2
        end
    end

    Y
end


# Solve for multiple RHS
# TODO: Generalize def
function cont_constrained_lyapuptriang{N,OSS<:OperatorSchur,PF<:ProductFun}(::Type{N},
                                                                            OS::PDEOperatorSchur{OSS},
                                                                            Gx,
                                                                            F::Array{PF,2};kwds...)
    n = min(size(OS.S.T,2),max(mapreduce(Fk->size(Fk,2),max,F),mapreduce(Fk->size(Fk,2),max,Gx)))
    rs=rangespace(OS)

    Fm=size(F,2)

    Y=Array(Fun{typeof(domainspace(OS,1)),N},n,Fm) # the solution
    PY=Array(Fun{typeof(rs[1]),N},n,Fm) # the first x-op times the solution
    SY=Array(Fun{typeof(rs[1]),N},n,Fm) # the second x-op times the solution


    k=n
    m=n  # max length



    nbcs=size(first(Gx),1)
    rhs=Array(Any,nbcs+1,Fm)

    TT=isempty(OS.Bx)?eltype(OS):promote_type(eltype(OS),mapreduce(eltype,promote_type,OS.Bx))




     while k≥1
        @assert k==1 || (OS.S.T[k,k-1] == 0 && OS.S.R[k,k-1] == 0        )

        for j=1:Fm
            rhs[1:nbcs,j]=Gx[j][:,k]
            rhs[end,j]=F[j].coefficients[k]
            if k < n
                for kk=k+1:n
                    axpy!(-OS.S.R[k,kk],PY[kk,j],rhs[end,j]) # equivalent to X+=a*Y
                    axpy!(-OS.S.T[k,kk],SY[kk,j],rhs[end,j])
                end
            end
        end
        ops=OS.Rdiags[k]
        Y[k,:]=vec(linsolve(ops,rhs;kwds...))
        for j=1:Fm
            chop!(Y[k,j],eps())
        end

        if k > 1
            PY[k,:]=vec(OS.Lx*Y[k:k,:])
            SY[k,:]=vec(OS.Mx*Y[k:k,:])  # L*Array returns an array space op
        end
        k-=1
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



cont_constrained_lyap{OSS<:OperatorSchur}(OS::PDEOperatorSchur{OSS},
                                          Gx::Vector,
                                          Gyin::Vector,
                                          F::Fun;kwds...)=cont_constrained_lyap(OS,Gx,Gyin,ProductFun(F);kwds...)
function cont_constrained_lyap{OSS<:OperatorSchur}(OS::PDEOperatorSchur{OSS},
                                                   Gx::Vector,
                                                   Gyin::Vector,
                                                   F::ProductFun;kwds...)
    Gy=regularize_bcs(OS.S,Gyin)
    F=pad!(F,:,min(size(OS.S.Q,1),size(F,2)))

    F=cont_reduce_dofs!(OS.S,OS.Lx,OS.Mx,Gy,F)

     # Q2 says how to rearrange the columns of F so that the operator is upper triangular
    Q2 = OS.S.Q
    F=ProductFun(Q2[1:length(F.coefficients),:]'*F.coefficients,space(F))


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

    Y=cont_constrained_lyapuptriang(OS,Gx,F;kwds...)
    # Y is a Vector{Fun}, so that Y[k][j] corresponds to matrix element M[k,j]
    # This means acting on Y is acting on *columns* of M

    X22=OS.S.Z*Y  #think of it as transpose

    X11=Gy-OS.S.bcs[:,Ky+1:end]*X22
    X=[X11;X22]
    X=OS.S.bcP*X        # this is equivalent to acting on columns by P'
end

function cont_constrained_lyap{OSS<:OperatorSchur,PF<:ProductFun}(OS::PDEOperatorSchur{OSS},
                                                                  Gx::Array,
                                                                  Gyin::Array,
                                                                  F::Array{PF};kwds...)
    @assert size(F,1)==1

    Gy=regularize_bcs(OS.S,Gyin)
    cont_reduce_dofs!(OS.S,OS.Lx,OS.Mx,Gy,F)
    Q2 = OS.S.Q
    for k=1:size(F,2)
        F[1,k]=ProductFun(Q2[1:length(F[1,k].coefficients),:].'*F[1,k].coefficients,space(F[1,k]))
    end
    ny=size(OS.S,2)
    Ky=numbcs(OS.S)

    if !isempty(Gx)
        GxM=Array(Matrix{Float64},size(Gx,2))
        for k=1:size(Gx,2)
            # bcP recombines boundary conditions
            GxM[k]=pad(coefficients(Gx[:,k]).',:,ny)*OS.S.bcP
            # remove unused DOFs and rearrange columns
            GxM[k]=GxM[k][:,Ky+1:end]*OS.S.Z
        end
    else
        GxM=[]
    end
    nx=10000
    Y=cont_constrained_lyapuptriang(Float64,OS,GxM,F;kwds...)
    #TODO: Ky or Kx?
    X=Array(eltype(Y),size(Y,1)+Ky,size(Y,2))
    for k=1:size(Y,2)
        X[Ky+1:end,k]=OS.S.Z*Y[:,k]
        X[1:Ky,k]=Gy[:,k]-OS.S.bcs[:,Ky+1:end]*X[Ky+1:end,k]
        X[:,k]=OS.S.bcP*X[:,k]
    end

    X
end


cont_constrained_lyap{OSS<:OperatorSchur,PF<:Fun}(OS::PDEOperatorSchur{OSS},
                                      Gx::Array,
                                      Gyin::Array,
                                      F::Array{PF};kwds...)=cont_constrained_lyap(OS,Gx,Gyin,map(ProductFun,F);kwds...)
