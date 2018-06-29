import Base.BLAS.@blasfunc
import LinearAlgebra: chkstride1, BlasInt
import LinearAlgebra.LAPACK.chklapackerror

export QuotientSpace

struct QuotientSpace{S,O<:Operator,DD,T,RT} <: Space{DD,T}
    space::S
    bcs::O
    A::Matrix{T}
    x::Vector{T}
    b::Vector{T}
    U::Matrix{T}
    Σ::Diagonal{T}
    VT::Matrix{T}
    work::Vector{T}
    iwork::Vector{BlasInt}
    rwork::Vector{RT}
    info::Ref{BlasInt}
    function QuotientSpace{S,O,DD,T,RT}(space::S, bcs::O) where {S,O,DD,T,RT}
        n = size(bcs, 1)
        @assert isfinite(n)
        A = zeros(T, n, n)
        x = zeros(T, n)
        b = zeros(T, n)
        U = zeros(T, n, n)
        Σ = Diagonal(zeros(T, n))
        VT = zeros(T, n, n)
        work = Vector{T}(max((3n+7)*n,68))
        iwork = Vector{BlasInt}(8*n)
        rwork = Vector{RT}((5*n+7)*n)
        info = Ref{BlasInt}()
        new{S,O,DD,T,RT}(space, bcs, A, x, b, U, Σ, VT, work, iwork, rwork, info)
    end
end

QuotientSpace(sp::Space, bcs::Operator) = QuotientSpace{typeof(sp), typeof(bcs), domaintype(sp), rangetype(sp), real(rangetype(sp))}(sp, bcs)

QuotientSpace(bcs::Operator) = QuotientSpace(domainspace(bcs), bcs)

domain(QS::QuotientSpace) = domain(QS.space)
dimension(QS::QuotientSpace) = dimension(QS.space)

Conversion(Q::QuotientSpace{SP}, S::SP) where {SP<:Space} = ConcreteConversion(Q, S)

bandinds(C::ConcreteConversion{QuotientSpace{SP,O,DD,T,RT},SP}) where {SP,O,DD,T,RT} = 0, size(C.domainspace.A, 1)

function getindex(C::ConcreteConversion{QuotientSpace{SP,O,DD,T,RT},SP}, i::Integer, j::Integer) where {SP,O,DD,T,RT}
    sp = domainspace(C)
    n = size(sp.A, 1)
    A = sp.A
    x = sp.x
    b = sp.b
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
    elseif i ≤ n && j ≤ n
        zero(T)
    elseif j-n ≤ i < j
        for jj = 1:n, ii = 1:n
            A[ii,jj] = B[ii,j+jj-n-1]
        end
        for ii = 1:n
            b[ii] = -B[ii,j]
        end
        in_place_gesdd!(A, U, Σ.diag, VT, work, iwork, rwork, info)
        mul!(x, transpose(U), b)
        mul!(b, pinv!(Σ), x)
        mul!(x, transpose(VT), b)
        x[i+n-j+1]
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
            job = 'A'
            lwork  = BlasInt(-1)
            cmplx  = $elty != $relty
            for i = 1:2
                if cmplx
                    ccall((@blasfunc($gesdd), liblapack), Nothing,
                          (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                           Ptr{BlasInt}, Ptr{$relty}, Ptr{$elty}, Ptr{BlasInt},
                           Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                           Ptr{$relty}, Ptr{BlasInt}, Ptr{BlasInt}),
                          &job, &m, &n, A, &max(1,stride(A,2)), S, U, &max(1,stride(U,2)), VT, &max(1,stride(VT,2)),
                          work, &lwork, rwork, iwork, info)
                else
                    ccall((@blasfunc($gesdd), liblapack), Nothing,
                          (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                           Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt},
                           Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                           Ptr{BlasInt}, Ptr{BlasInt}),
                          &job, &m, &n, A, &max(1,stride(A,2)), S, U, &max(1,stride(U,2)), VT, &max(1,stride(VT,2)),
                          work, &lwork, iwork, info)
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
