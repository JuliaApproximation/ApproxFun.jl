
##
# These routines simplify the multiplication of two operators
# they are currently used for dekron of PDEs to make sure equivalent operators are
# recognized as the same
# I put this in Extras as it's currently a "hack"
#
# In the future, product rule/etc. could be implemented, and some sort of ordering
# to recognize M[f,C^(1)] Conversion(T,C^(1))==Conversion(T,C^(1))M[f,T]
##

simplifytimes(A::Operator,B::Operator) = isconstop(A)||isconstop(B) ? A*B : [A;B]
simplifytimes(A::Operator,B::AbstractVector) = [simplifytimes(A,B[1]);B[2:end]]
simplifytimes(A::AbstractVector,B::Operator) = [A[1:end-1];simplifytimes(A[end],B)]
function simplifytimes(A::ConstantOperator,B::ConstantOperator)
    @assert A.space == B.space
    ConstantOperator(Number(A)*Number(B),A.space)
end
function simplifytimes(A::SpaceOperator,B::SpaceOperator)
   @assert domainspace(A)==rangespace(B)
    SpaceOperator(simplifytimes(A.op,B.op),domainspace(B),rangespace(A))
end

simplifytimes(A::Conversion,B::Conversion) = Conversion(domainspace(B),rangespace(A))

simplify(A)=A

function simplify(A::TimesOperator)
    ret=foldr(simplifytimes,A.ops)
    if isa(ret,Vector)
        @assert length(ret)>1
        TimesOperator(ret)
    else
        ret
    end
end
