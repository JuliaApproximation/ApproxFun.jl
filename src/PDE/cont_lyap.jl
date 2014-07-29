function adaptiveplus(f::Vector,g::Vector)
    if length(f)>length(g)
        ret=copy(f)
        ret[1:length(g)]+=g
    else
        ret=copy(g)
        ret[1:length(f)]+=f
    end
    
    ret
end
        
adaptiveplus(f,g,h)=adaptiveplus(adaptiveplus(f,g),h)


#Use XR' = G' = [G1 G2 G3...] to reduce columns of A in
# MXA' + *X* =F
# here G is a vector of IFuns


function cont_reduce_dofs{T<:IFun}( R,G::Vector{T}, A::Array, M::Operator, F::Array )
    if length(R) > 0
        # first multiply to get MXR' = M*G' = [M*G1 M*G2 ...]
        # then kill the row by subtracting
        # MXR'[:,k]*A'[k,:]  from MXA'
        # i.e., subtacting A[:,k]*R[k,:] from A
        # and M*G'[:,k]*A'[k,:] from F
        # i.e. M*G[k]*A[:,k]' from 
        
        for k = 1:size(R,1)
            MG = M*G[k].coefficients         # coefficients in the range space of M      
            F = F - pad(MG*A[:,k].',size(F,1),size(F,2))
            A = A - A[:,k]*R[k,:]
        end
    end
        
    A, F
end



# Solve Bx*Y=Gx and P*Y*R' + S*Y*T' = F 
# where R and T are upper triangular
function cont_constrained_lyapuptriang{N}(Bx,Gx,P,R,S,T,F::Array{N})
    n = size(T,2)
    Y=Array(Vector{N},n)
    PY=Array(Vector{N},n)
    SY=Array(Vector{N},n)

    k=n
    while k>1
        if T[k,k-1] == 0 && R[k,k-1] == 0        
            rhs = F[:,k]
            if k < n
                for j=k+1:n
                    rhs= adaptiveplus(rhs,-R[k,j]*PY[j],-T[k,j]*SY[j])
                end
            end

            op=(R[k,k]*P + T[k,k]*S);
            Y[k]=chop!([Bx,op]\[Gx[:,k],rhs],eps());
            
            PY[k]=P*Y[k];
            SY[k]=S*Y[k];   
            k-=1
        else
            error("non-upper triangular not implemented")
        end
    end

    if k == 1
        rhs = F[:,k]

        for j=2:n
            rhs= adaptiveplus(rhs,-R[k,j]*PY[j],-T[k,j]*SY[j])
        end

        Y[k]=chop!([Bx,(R[k,k]*P + T[k,k]*S)]\[Gx[:,k],rhs],eps());
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


function cont_constrained_lyap(Bxin,Byin,Lin,Min,F,ny)
    Xop=promotespaces([Lin[1],Min[1]])
    Lx=SavedBandedOperator(Xop[1]);Mx=SavedBandedOperator(Xop[2])

    #discretize in Y
    By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny) 
    

    Ry,Gy,Ly,My,Py=regularize_bcs(By,Gy,Ly,My) 
    Ly,F = cont_reduce_dofs(Ry,Gy,Ly,Lx,F)     

    Ky = size(By,1)
    B=Ly[:,Ky+1:end]
    D=My[:,Ky+1:end]
    BD=schurfact(full(B),full(D))
    Q2=BD[:left];Z2=BD[:right]
    R=BD[:S]; T=BD[:T]

    F=pad(F,size(F,1),size(Q2,1))*Q2
    

    ## we've discretized, in y, and rhs for Bx is a function of y
    # so we need to discetize it as well
    Gx=toarray(Bxin[2],ny)
    # remove unused DOFs and rearrange columns
    Gx=Gx[:,Ky+1:end]*Z2


    Y=cont_constrained_lyapuptriang(Bxin[1],Gx,Lx,R,Mx,T,F)
    m=mapreduce(length,max,Y)
    Y22=hcat([pad(Y[k],m) for k=1:length(Y)]...)
    X22=Y22*Z2'
    X11=pad(Gy',size(X22,1),2)-X22*Ry[:,Ky+1:end]'
    [X11 X22]
end


