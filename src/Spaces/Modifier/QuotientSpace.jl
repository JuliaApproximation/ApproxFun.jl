import LinearAlgebra: LU, checknonsingular

export QuotientSpace

struct QuotientSpace{S,O,D,R} <: Space{D,R}
    space::S
    bcs::O
    F::LU{R,Matrix{R}}
    x::Vector{R}
    function QuotientSpace{S,O,D,R}(space::S, bcs::O) where {S,O,D,R}
        n = size(bcs, 1)
        @assert isfinite(n)
        new{S,O,D,R}(space, bcs, lu!(Matrix{R}(I, n, n)), zeros(R, n))
    end
end

QuotientSpace(sp::Space{D,R}, bcs::Operator) where {D,R} = QuotientSpace{typeof(sp), typeof(bcs), D, R}(sp, bcs)
QuotientSpace(bcs::Operator) = QuotientSpace(domainspace(bcs), bcs)

domain(QS::QuotientSpace) = domain(QS.space)
canonicalspace(QS::QuotientSpace) = QS.space

spacescompatible(a::QuotientSpace,b::QuotientSpace) = spacescompatible(a.space,b.space)
hasconversion(a::Space,b::QuotientSpace) = hasconversion(a,b.space)
hasconversion(a::QuotientSpace,b::Space) = hasconversion(a.space,b)
hasconversion(a::QuotientSpace,b::QuotientSpace) = hasconversion(a.space,b)


Conversion(Q::QuotientSpace{S}, sp::S) where S<:Space = ConcreteConversion(Q, sp)
@inline bandwidths(C::ConcreteConversion{QuotientSpace{S,O,D,R},S}) where {S,O,D,R} = (size(C.domainspace.F, 1), 0)

function getindex(C::ConcreteConversion{QuotientSpace{S,O,D,R},S}, i::Integer, j::Integer) where {S,O,D,R}
    sp = domainspace(C)
    B = sp.bcs
    F = sp.F
    A = F.factors
    n = size(A, 1)
    x = sp.x
    if i == j
        return one(R)
    elseif j < i ≤ j+n
        @inbounds for jj = 1:n, ii = 1:n
            A[ii,jj] = B[ii,j+jj]
        end
        @inbounds for ii = 1:n
            x[ii] = -B[ii,j]
        end
        if norm(x) > 8*norm(A)*eps(R)
            mutable_lu!(F)
            ldiv!(F, x)
        end
        return @inbounds x[i-j]
    else
        return zero(R)
    end
end

# This function is modified from Julia's lu.jl in stdlib to work on a single `LU`.
# License is MIT: https://julialang.org/license
function mutable_lu!(F::LU{T}, ::Val{Pivot} = Val(true);
                         check::Bool = true) where {T,Pivot}
    A = F.factors
    m, n = size(A)
    minmn = min(m,n)
    info = 0
    ipiv = F.ipiv
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if Pivot
                amax = abs(zero(T))
                for i = k:m
                    absi = abs(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    check && checknonsingular(info)
    return F
end



import LinearAlgebra.BLAS.@blasfunc
import LinearAlgebra: chkstride1, BlasInt
import LinearAlgebra.LAPACK.chklapackerror

export PathologicalQuotientSpace

struct PathologicalQuotientSpace{S,O<:Operator,DD,T,RT} <: Space{DD,T}
    space::S
    bcs::O
    A::Matrix{T}
    x::Vector{T}
    b::Vector{T}
    c::Vector{T}
    U::Matrix{T}
    Σ::Diagonal{T,Vector{T}}
    VT::Matrix{T}
    work::Vector{T}
    iwork::Vector{BlasInt}
    rwork::Vector{RT}
    info::Ref{BlasInt}
    function PathologicalQuotientSpace{S,O,DD,T,RT}(space::S, bcs::O) where {S,O,DD,T,RT}
        n = size(bcs, 1)
        @assert isfinite(n)
        A = zeros(T, n, 2n)
        x = zeros(T, 2n)
        b = zeros(T, n)
        c = zeros(T, n)
        U = zeros(T, n, n)
        Σ = Diagonal(zeros(T, n))
        VT = zeros(T, n, 2n)
        work = Vector{T}(undef,max((6n+7)*2n,68))
        iwork = Vector{BlasInt}(undef,16n)
        rwork = Vector{RT}(undef,(10n+7)*2n)
        info = Ref{BlasInt}()
        new{S,O,DD,T,RT}(space, bcs, A, x, b, c, U, Σ, VT, work, iwork, rwork, info)
    end
end

PathologicalQuotientSpace(sp::Space, bcs::Operator) = PathologicalQuotientSpace{typeof(sp), typeof(bcs), domaintype(sp), rangetype(sp), real(rangetype(sp))}(sp, bcs)

PathologicalQuotientSpace(bcs::Operator) = PathologicalQuotientSpace(domainspace(bcs), bcs)

domain(QS::PathologicalQuotientSpace) = domain(QS.space)
dimension(QS::PathologicalQuotientSpace) = dimension(QS.space)

Conversion(Q::PathologicalQuotientSpace{SP}, S::SP) where {SP<:Space} = ConcreteConversion(Q, S)

bandwidths(C::ConcreteConversion{PathologicalQuotientSpace{SP,O,DD,T,RT},SP}) where {SP,O,DD,T,RT} = (2size(C.domainspace.A, 1), 0)

function getindex(C::ConcreteConversion{PathologicalQuotientSpace{SP,O,DD,T,RT},SP}, i::Integer, j::Integer) where {SP,O,DD,T,RT}
    sp = domainspace(C)
    n = size(sp.A, 1)
    A = sp.A
    x = sp.x
    b = sp.b
    c = sp.c
    B = sp.bcs
    U = sp.U
    Σ = sp.Σ
    VT = sp.VT
    work = sp.work
    iwork = sp.iwork
    rwork = sp.rwork
    info = sp.info

    if i == j
        one(T)
    elseif j < i ≤ j+2n
        for jj = 1:2n, ii = 1:n
            A[ii,jj] = B[ii,j+jj]
        end
        for ii = 1:n
            b[ii] = -B[ii,j]
        end
        in_place_gesdd!(A, U, Σ.diag, VT, work, iwork, rwork, info)
        mul!(c, transpose(U), b)
        mul!(b, pinv!(Σ), c)
        mul!(x, transpose(VT), b)
        x[i-j]
    else
        zero(T)
    end
end

function pinv!(D::Diagonal{T}) where T
    d = D.diag
    @inbounds @simd for i = 1:length(d)
        isfinite(inv(d[i])) ? d[i]=inv(d[i]) : d[i]=zero(T)
    end
    D
end

for (gesdd, elty, relty) in ((:dgesdd_,:Float64,:Float64),
                      (:sgesdd_,:Float32,:Float32),
                      (:zgesdd_,:ComplexF64,:Float64),
                      (:cgesdd_,:ComplexF32,:Float32))
    @eval begin
        #    SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
        #                   LWORK, IWORK, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          JOBZ
        #      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
        #*     ..
        #*     .. Array Arguments ..
        #      INTEGER            IWORK( * )
        #      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
        #                        VT( LDVT, * ), WORK( * )
        function in_place_gesdd!(A::Matrix{$elty}, U::Matrix{$elty}, S::Vector{$elty}, VT::Matrix{$elty}, work::Vector{$elty}, iwork::Vector{BlasInt}, rwork::Vector{$relty}, info::Ref{BlasInt})
            chkstride1(A)
            m, n   = size(A)
            minmn  = min(m, n)
            job = 'S'
            lwork  = BlasInt(-1)
            cmplx  = $elty != $relty
            for i = 1:2
                if cmplx
                    ccall((@blasfunc($gesdd), liblapack), Nothing,
                          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$relty}, Ptr{BlasInt}, Ptr{BlasInt}),
                          job, m, n, A,
                          max(1,stride(A,2)), S, U, max(1,stride(U,2)),
                          VT, max(1,stride(VT,2)), work, lwork,
                          rwork, iwork, info)
                else
                    ccall((@blasfunc($gesdd), liblapack), Nothing,
                          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{BlasInt}, Ptr{BlasInt}),
                          job, m, n, A,
                          max(1,stride(A,2)), S, U, max(1,stride(U,2)),
                          VT, max(1,stride(VT,2)), work, lwork,
                          iwork, info)
                end
                chklapackerror(info[])
                if i == 1
                    # Work around issue with truncated Float32 representation of lwork in
                    # sgesdd by using nextfloat. See
                    # http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=4587&p=11036&hilit=sgesdd#p11036
                    # and
                    # https://github.com/scipy/scipy/issues/5401
                    lwork = round(BlasInt, nextfloat(real(work[1])))
                    #work = Vector{$elty}(lwork)
                end
            end
        end
    end
end
