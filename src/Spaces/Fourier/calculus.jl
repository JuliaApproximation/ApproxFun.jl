##Differentiation and integration


#differentiate(f::Fun{Laurent})=Fun(interlace(fourierdiff(domain(f),deinterlace(f.coefficients))),f.space)
Base.sum(f::Fun{Laurent})=fouriersum(domain(f),deinterlace(f.coefficients))
function integrate(f::Fun{Hardy{false}})
    if isa(domain(f),Circle) # drop -1 term if zero and try again
        integrate(Fun(f,DropSpace(space(f),1)))
    else  # Probably periodic itnerval
        Integral(space(f))*f
    end
end

function integrate(f::Fun{Taylor})
    if isa(domain(f),Circle)
        Integral(space(f))*f
    else  # Probably periodic itnerval  drop constant term if zero
        integrate(Fun(f,DropSpace(space(f),1)))
    end
end


function integrate(f::Fun{CosSpace})
    if isa(domain(f),Circle)
        error("Integrate not implemented for CosSpace on Circle")
    else  # Probably periodic itnerval, drop constant term if zero
        integrate(Fun(f,DropSpace(space(f),1)))
    end
end

function integrate(f::Fun{SinSpace})
    if isa(domain(f),Circle) # drop term containing z^(-1)
        integrate(Fun(f,DropSpace(space(f),1)))
    else  # Probably periodic itnerval\
        Integral(space(f))*f
    end
end

#TODO: Hack to make sure Fourier maps to Fourier
for OP in (:differentiate,:integrate)
    @eval $OP{T}(f::Fun{Fourier,T})=$OP(vec(f,2))⊕$OP(vec(f,1))
end
# 
# 
# fourierdiff(d::PeriodicInterval,cfs::ShiftVector)=tocanonicalD(d,0)*ShiftVector(1.im*[firstindex(cfs):-1],1.im*[0:lastindex(cfs)]).*cfs
# 
# 
# 
# 
# function fourierintegrate(d::PeriodicInterval,cfs::ShiftVector)
#     tol = 10eps()
#     @assert abs(cfs[0]) < tol
#     
#     ##TODO: mapped domains
#     
#     @assert d.a ==-π
#     @assert d.b ==π        
#     ShiftVector(-1im*cfs[firstindex(cfs):-1]./[firstindex(cfs):-1],
#                 [0,(-1im*cfs[1:lastindex(cfs)]./[1:lastindex(cfs)])])
# end
# 

# 
# 
# 
# function fourierdiff(d::Circle,cfs::ShiftVector)
#         ##TODO: general radii
#         @assert d.radius == 1.
#         @assert d.center == 0
# 
#         # Now shift everything by one
#         ShiftVector(
#                         [cfs[firstindex(cfs):-1].*[firstindex(cfs):-1],0],
#                         cfs[1:lastindex(cfs)].*[1:lastindex(cfs)]
#                         )
# end
# 
# 
# 
# function fourierintegrate(d::Circle,cfs::ShiftVector)
#     tol = 10eps()
#     @assert abs(cfs[-1]) < tol        
#     ##TODO: general radii        
#     @assert d.radius == 1.
#     @assert d.center == 0        
#     
#     # Now shift everything by one
#     ShiftVector(
#                     [cfs[firstindex(cfs):-1]./[firstindex(cfs):-1]],
#                     [0,(cfs[0:lastindex(cfs)]./[1:lastindex(cfs)+1])]
#                     )
# end
# 


fouriersum(d::PeriodicInterval,cfs::ShiftVector)=cfs[0].*length(d)

function fouriersum{T}(d::Circle,cfs::ShiftVector{T})
    @assert d.radius == 1.
    @assert d.center == 0   
    if firstindex(cfs) <= -1
        cfs[-1]
    else
        zero(T)
    end
end

