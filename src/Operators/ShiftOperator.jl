## BandedShiftOperator allows automatic interlacing

abstract BandedShiftOperator{T} <: BandedOperator{T}
abstract ShiftFunctional{T} <: Functional{T}


function bandinds(b::BandedShiftOperator)
    bi=shiftbandinds(b)
    m=max(-bi[1],bi[2])
    -2m,2m
end

shiftbandrange(b::BandedShiftOperator)=Range1(shiftbandinds(b)...)



shiftShiftArray{T<:Number}(B::BandedShiftOperator{T},k::Range1,j::Range1)=shiftaddentries!(B,sazeros(T,k,j),k)
shiftShiftArray(B::BandedShiftOperator,k::Range1)=shiftShiftArray(B,k,shiftbandrange(B))
shiftBandedArray(B::BandedShiftOperator,k::Range1)=shiftBandedArray(B,k,(k[1]+shiftbandinds(B)[1]):(k[end]+shiftbandinds(B)[end]))
shiftBandedArray(B::BandedShiftOperator,k::Range1,cs)=BandedArray(shiftShiftArray(B,k,shiftbandrange(B)),cs)


# BandedShiftOperator overrides shiftaddentries!


shiftfirstrw(rs,ri,k::Integer)=rs>0?fld(k-ri+rs-1,rs):fld(k-ri,rs)
shiftfirstrw(S,k::Integer)=firstrw(S.rowstride,S.rowindex,k)

#Last index below
shiftlastrw(rs,ri,k::Integer)=rs>0?fld(k-ri,rs):fld(k-ri+rs+1,rs)

function shiftdivrowrange(rs,ri,r)
    if rs > 0
        shiftfirstrw(rs,ri,r[1]):shiftlastrw(rs,ri,r[end])
    else #neg neg
        shiftlastrw(rs,ri,r[end]):shiftfirstrw(rs,ri,r[1])
    end
end


function shift_stride_pospos_addentries!(ri,ci,rs,cs,S,A,kr::Range)
    r1=shiftdivrowrange(rs,ri,kr)

    B1=shiftBandedArray(S,r1)
    B=BandedArray(A)
    
    for k=r1, j=columnrange(B1.data)+k
        B[rs*k + ri,cs*j + ci] += B1.data[k,j-k]
    end
    
    A
end

function shift_stride_posneg_addentries!(ri,ci,rs,cs,S,A,kr::Range)
    r1=shiftdivrowrange(rs,ri,kr)
    B1=shiftShiftArray(S,r1)
    B=BandedArray(A)
    
    br=shiftbandrange(S)
    
    for k=r1, j=br
        if cs*(j+k) + ci > 0 && rs*k + ri > 0
            B[rs*k + ri,cs*(j+k) + ci] += B1[k,j]
        end
    end

    
    A
end

function shift_stride_addentries!(ri,ci,rs,cs,S,A,kr)
    if rs > 0 && cs > 0
        shift_stride_pospos_addentries!(ri,ci,rs,cs,S,A,kr)
    elseif rs > 0
        shift_stride_posneg_addentries!(ri,ci,rs,cs,S,A,kr)    
    elseif cs > 0
        shift_stride_posneg_addentries!(ri,ci,rs,cs,S,A,kr)            
    else #neg neg
        shift_stride_pospos_addentries!(ri,ci,rs,cs,S,A,kr)            
    end
end


function addentries!(L::BandedShiftOperator,A,kr::Range)
    shift_stride_addentries!(1,1,2,2,L,A,kr)
    shift_stride_addentries!(0,0,-2,-2,L,A,kr)
    shift_stride_addentries!(0,1,-2,2,L,A,kr)
    shift_stride_addentries!(1,0,2,-2,L,A,kr)    
end

