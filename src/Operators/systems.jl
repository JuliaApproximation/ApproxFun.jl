##Operators

for op in (:DerivativeOperator,:IntegrationOperator)
    @eval begin
        function ($op){T<:IntervalDomain}(d::Vector{T})
            n=length(d)
            R=zeros(Operator,n,n)
            for k=1:n
                R[k,k]=$op(d[k])
            end
            
            R
        end
    end
end




        function EvaluationFunctional{T<:IntervalDomain}(d::Vector{T},x...)
            n=length(d)
            R=zeros(Operator,n,n)
            for k=1:n
                R[k,k]=EvaluationFunctional(d[k],x...)
            end
            
            R
        end



