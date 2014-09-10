## Multiplication of operator * fun

# function ultraiconversion{T}(g::Vector{T},m::Integer)
#     if m==0 
#         g 
#     elseif m==1
#         ultraiconversion(g)
#     else        # here domain is arbitrary
#         backsubstitution!(MutableAlmostBandedOperator(Operator[USConversionOperator(0:m,Interval())]),copy(g))::Vector{T}
#     end
# end
# function ultraconversion{T}(g::Vector{T},m::Integer)
#     if m==0
#         g 
#     elseif m==1
#         ultraconversion(g)
#     else # here domain is arbitrary
#         (USConversionOperator(0:m,Interval())*g)::Vector{T}
#     end
# end


function maxspace{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder > border?a:b
end

function minspace{aorder,border}(a::UltrasphericalSpace{aorder},b::UltrasphericalSpace{border})
    @assert domainscompatible(a,b)
    
    aorder < border?a:b
end


spaceconversion(g::Vector,::UltrasphericalSpace{1},::ChebyshevSpace)=ultraiconversion(g)
spaceconversion(g::Vector,::ChebyshevSpace,::UltrasphericalSpace{1})=ultraconversion(g)

