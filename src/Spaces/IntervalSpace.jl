export continuity



Space(d::IntervalDomain) = Chebyshev(d)

identity_fun(d::Interval) = Fun(eltype(d)[(d.b+d.a)/2,(d.b-d.a)/2],Chebyshev(d))


## Calculus



# the default domain space is higher to avoid negative ultraspherical spaces
Integral(d::IntervalDomain,n::Integer) = Integral(Ultraspherical(1,d),n)

for Func in (:DefiniteIntegral,:DefiniteLineIntegral)
    @eval begin
        #TODO: this may be misleading
        $Func(d::IntervalDomain) = $Func(JacobiWeight(-.5,-.5,Chebyshev(d)))
        function $Func(α::Number,β::Number,d::IntervalDomain)
            @assert α == β
            @assert round(Int,α+.5) == α+.5
            @assert round(Int,α+.5) >= 0
            $Func(JacobiWeight(α,β,Ultraspherical(round(Int,α+.5),d)))
        end
        $Func(α::Number,β::Number) = $Func(α,β,Interval())
    end
end



## Evaluation








function continuity(d::Union{Vector,Tuple},order::Integer)

    m=length(d)
    B=zeros(Operator{mapreduce(eltype,promote_type,d)},m-1,m)

    for k=1:m-1
        B[k,k]=Evaluation(d[k],true,order)
        B[k,k+1]=-Evaluation(d[k+1],false,order)
    end
    B
end

function continuity(d::Union{Vector,Tuple},kr::UnitRange)
    @assert first(kr)==0
    m=length(d)
    B=zeros(Operator{mapreduce(eltype,promote_type,d)},length(kr)*(m-1),m)
    for r in kr
        B[(m-1)*r+1:(m-1)*(r+1),:]=continuity(d,r)
    end
    B
end


for OP in (:dirichlet,:neumann,:periodic)
    @eval $OP(d::Tuple)=$OP([d...])
end

for DT in (:IntervalDomain,:Space)
    @eval begin
        function dirichlet{T<:$DT}(d::Vector{T})
            m=length(d)
            B=zeros(Operator{mapreduce(eltype,promote_type,d)},2,m)
            B[1,1]=ldirichlet(d[1]);B[2,end]=rdirichlet(d[end])
            [B;
            continuity(d,0:1)]
        end

        function neumann{T<:$DT}(d::Vector{T})
            m=length(d)
            B=zeros(Operator{mapreduce(eltype,promote_type,d)},2,m)
            B[1,1]=ldirichlet(d[1]);B[2,end]=rdirichlet(d[end])
            [B;
            continuity(d,0:1)]
        end


        function periodic{T<:$DT}(d::Vector{T})
            m=length(d)
            B=zeros(Operator{mapreduce(eltype,promote_type,d)},2,m)
            B[1,1]=ldirichlet(d[1]);B[1,end]=-rdirichlet(d[end])
            B[2,1]=lneumann(d[1]);B[2,end]=-rneumann(d[end])
            [B;
            continuity(d,0:1)]
        end
    end
end
