

immutable Hermite{T} <: PolynomialSpace
    L::T
end

    spacescompatible(::Hermite,::Hermite)=true #TODO:L

    domain(::Hermite)=Line()
    Hermite()=Hermite(1.0)



    #####
    # recα/β/γ are given by
    #       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
    #####


    recα(::Type,::Hermite,k)=0;recβ(::Type,::Hermite,k)=0.5;recγ(::Type,::Hermite,k)=k-1

    bandinds{H<:Hermite}(D::Derivative{H})=0,D.order
    rangespace{H<:Hermite}(D::Derivative{H})=domainspace(D)
    function addentries!{H<:Hermite}(D::Derivative{H},A,kr::Range)
        if D.order==1
            for k=kr
               A[k,k+D.order]+=2k
            end
        elseif D.order==2
            for k=kr
               A[k,k+D.order]+=4k*(k+1)
            end
        else
            error("Not implemented")
        end
        A
    end


    identity_fun(sp::Hermite)=Fun([0.,0.5],sp)


# exp(-Lx^2)
immutable GaussWeight{S,T} <: WeightSpace
    space::S
    L::T
end
spacescompatible(a::GaussWeight,b::GaussWeight)=spacescompatible(a.space,b.space)&&isapprox(a.L,b.L)

function Derivative(sp::GaussWeight,k)
   if k==1
        x=Multiplication(Fun(identity,sp.space),sp.space)
        D=Derivative(sp.space)
        D2=D-(2sp.L)x
        DerivativeWrapper(SpaceOperator(D2,sp,GaussWeight(rangespace(D2),sp.L)),1)
    else
        D=Derivative(sp)
        DerivativeWrapper(TimesOperator(Derivative(rangespace(D),k-1).op,D.op),k)
    end
end

weight(H::GaussWeight,x)=exp(-H.L*x.^2)

Derivative(GaussWeight(Hermite(),0.5),1)[1:10,1:10]

