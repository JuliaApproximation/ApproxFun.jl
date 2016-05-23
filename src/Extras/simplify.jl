
##
# These routines simplify the multiplication of two operators
# they are currently used for dekron of PDEs to make sure equivalent operators are
# recognized as the same
# I put this in Extras as it's currently a "hack"
#
# In the future, product rule/etc. could be implemented, and some sort of ordering
# to recognize M[f,C^(1)] Conversion(T,C^(1))==Conversion(T,C^(1))M[f,T]
##

simplifytimes(A::BandedOperator,B::BandedOperator)=isconstop(A)||isconstop(B)?A*B:[A;B]
simplifytimes(A::BandedOperator,B::Vector)=[simplifytimes(A,B[1]);B[2:end]]
simplifytimes(A::Vector,B::BandedOperator)=[A[1:end-1];simplifytimes(A[end],B)]
simplifytimes(A::ConstantOperator,B::ConstantOperator)=ConstantOperator(A.c*B.c)
function simplifytimes(A::SpaceOperator,B::SpaceOperator)
   @assert domainspace(A)==rangespace(B)
    SpaceOperator(simplifytimes(A.op,B.op),domainspace(B),rangespace(A))
end

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
