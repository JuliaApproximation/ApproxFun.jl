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



function cont_reduce_dofs( R,G, Mx, My, F )
    if length(R) > 0
        GM = [My*(G.')[:,k] for k=1:size(G,1)]
        for k = 1:size(R,1)
            F = F - pad(Mx[k]*GM[k].',size(F,1),size(F,2))
                Mx = Mx - Mx[:,k]*R[k,:]
        end
    end
        
    Mx, F
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

function cont_constrained_lyap(Bxin,Byin,Lin,Min,F,ny)
    Xop=promotespaces([Lin[1],Min[1]])
    Lx=SavedBandedOperator(Xop[1]);Mx=SavedBandedOperator(Xop[2])
    Gx=Bxin[2]
    Bx=Bxin[1]
    
    By,Gy,Ly,My=pdetoarray(Byin,Lin[2],Min[2],ny)
    Ry,Gy,Ly,My,Py=regularize_bcs(By,Gy,Ly,My)
    Ly,F = cont_reduce_dofs(Ry,Gy,Ly,Lx,F.');    F = F.';

    Ky = size(By,1);
    B=Ly[:,Ky+1:end];
    D=My[:,Ky+1:end];
    BD=schurfact(full(B),full(D));
    Q2=BD[:left];Z2=BD[:right];
    R=BD[:S]; T=BD[:T];

    F=pad(F,size(F,1),size(Q2,1))*Q2
    
    Gx=pad(Gx,size(Gx,1),size(Z2,1)+2)[:,Ky+1:end]*Z2

    Y=cont_constrained_lyapuptriang(Bx,Gx,Lx,R,Mx,T,F);
    m=mapreduce(length,max,Y)
    Y22=hcat([pad(Y[k],m) for k=1:length(Y)]...);
    X22=Y22*Z2';
    X11=pad(Gy',size(X22,1),2)-X22*Ry[:,Ky+1:end]';
    [X11 X22]
end


