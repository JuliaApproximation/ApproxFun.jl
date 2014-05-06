

export PlusOperator,TimesOperator




type PlusOperator{T<:Number,B<:Operator} <: BandedOperator{T} 
    ops::Vector{B}
    
    PlusOperator{B}(ops::Vector{B})=new(promotespaces(ops))
end


typeofdata{T<:Number}(::Operator{T})=T
function PlusOperator{B<:Operator}(ops::Vector{B})
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    PlusOperator{T,B}(ops)
end


function domainspace(P::PlusOperator)
    for op in P.ops
        sp = domainspace(op)
        
        if sp != Any
            return sp
        end
    end
    
    Any
end

function rangespace(P::PlusOperator)
    for op in P.ops
        sp = rangespace(op)
        
        if sp != Any
            return sp
        end
    end
    
    Any
end

domain(P::PlusOperator)=commondomain(P.ops)

bandrange(P::PlusOperator)=mapreduce(op->bandrange(op)[1],min,P.ops):mapreduce(op->bandrange(op)[end],max,P.ops)


function addentries!(P::PlusOperator,A::ShiftArray,kr::Range1)
    for op in P.ops
        addentries!(op,A,kr)
    end
    
    A
end

+(A::PlusOperator,B::PlusOperator)=PlusOperator([A.ops,B.ops])
+(A::PlusOperator,B::Operator)=PlusOperator([A.ops,B])
+(A::Operator,B::PlusOperator)=PlusOperator([A,B.ops])
+(A::Operator,B::Operator)=PlusOperator([A,B])
+(A::Operator,f::AbstractFun)=A+MultiplicationOperator(f)
+(f::AbstractFun,A::Operator)=MultiplicationOperator(f)+A
-(A::Operator,f::AbstractFun)=A+MultiplicationOperator(-f)
-(f::AbstractFun,A::Operator)=MultiplicationOperator(f)-A
+(c::Number,A::Operator)=MultiplicationOperator(Fun(c,domain(A)))+A
.+(c::Number,A::Operator)=MultiplicationOperator(Fun(c,domain(A)))+A
+(A::Operator,c::Number)=A+MultiplicationOperator(Fun(c,domain(A)))
.+(A::Operator,c::Number)=A+MultiplicationOperator(Fun(c,domain(A)))
-(c::Number,A::Operator)=MultiplicationOperator(Fun(c,domain(A)))-A
-(A::Operator,c::Number)=A-MultiplicationOperator(Fun(c,domain(A)))



## Times Operator


type TimesOperator{T<:Number,B<:Operator} <: BandedOperator{T}
    ops::Vector{B}
    
    ##TODO: Probably should go other way
    function TimesOperator{B}(ops::Vector{B})
        for k=1:length(ops)-1
            ops[k]=promotedomainspace(ops[k],rangespace(ops[k+1]))
        end
        
        new(ops)
    end    
end

function TimesOperator{B<:Operator}(ops::Vector{B})
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    TimesOperator{T,B}(ops)
end


function domainspace(P::TimesOperator)
    for k=length(P.ops):-1:1
        sp = domainspace(P.ops[k])
        
        if sp != Any
            return sp
        end
    end
    
    Any
end

function rangespace(P::TimesOperator)
    for op in P.ops
        sp = rangespace(op)
        
        if sp != Any
            return sp
        end
    end
    
    Any
end

domain(P::TimesOperator)=commondomain(P.ops)


bandrange(P::TimesOperator)=mapreduce(x->bandrange(x)[1],+,P.ops):mapreduce(x->bandrange(x)[end],+,P.ops)


##TODO: We keep this around because its faster
## need to unify 
function old_addentries!{T<:Number,B}(P::TimesOperator{T,B},A::ShiftArray,kr::Range1)
    cr = columnrange(A)
    br = bandrange(P)

    kre=kr[1]:(kr[end]+mapreduce(x->bandrange(x)[end],+,P.ops[1:end-1]))

    Z = ShiftArray(P.ops[end],kre,br)
    
    for j=length(P.ops)-1:-1:2
        krr=kr[1]:(kr[end]+mapreduce(x->bandrange(x)[end],+,P.ops[1:j-1]))    
        multiplyentries!(P.ops[j],Z,krr)
    end
    
    multiplyentries!(P.ops[1],Z,kr)    
    
    for k=kr,j=max(cr[1],br[1]):min(cr[end],br[end])
        A[k,j] += Z[k,j]
    end
    
    A
end

function new_addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    krl=Array(Range1,length(P.ops))
    
    krl[1]=kr
    
    for m=1:length(P.ops)-1
        br=bandrange(P.ops[m])
         krl[m+1]=(br[1] + krl[m][1]):(br[end] + krl[m][end])
    end
    
    BA=BandedArray(P.ops[end],krl[end])
    
    for m=(length(P.ops)-1):-1:1
      BA=BandedArray(P.ops[m],krl[m],krl[m+1])*BA
    end
    
    for k=kr,j=columnrange(BA.data)
        A[k,j] += BA.data[k,j]
    end
    
    A
end

function addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    if all(f->typeof(f)<:ConversionOperator,P.ops[1:end-1])
        old_addentries!(P,A,kr)
    else
        new_addentries!(P,A,kr)
    end
end

*(A::TimesOperator,B::TimesOperator)=TimesOperator([A.ops,B.ops])
*(A::TimesOperator,B::Operator)=TimesOperator([A.ops,B])
*(A::Operator,B::TimesOperator)=TimesOperator([A,B.ops])
*(A::Operator,B::Operator)=TimesOperator([A,B])


-(A::Operator)=ConstantOperator(-1.)*A
-(A::Operator,B::Operator)=A+(-B)

*(f::IFun,A::BandedOperator)=MultiplicationOperator(f,rangespace(A).order)*A
*(f::FFun,A::Operator)=MultiplicationOperator(f)*A
*(c::Number,A::Operator)=ConstantOperator(c)*A



