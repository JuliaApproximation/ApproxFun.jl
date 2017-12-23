import Base.BLAS: @blasfunc, libblas
import Base.LAPACK: liblapack
# Level 2
## mv
### gemv

gemv!(x...) = BLAS.gemv!(x...)

for (fname, elty) in ((:dgemv_,:Float64),
                      (:sgemv_,:Float32),
                      (:zgemv_,:Complex128),
                      (:cgemv_,:Complex64))
    @eval begin
             #SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
             #*     .. Scalar Arguments ..
             #      DOUBLE PRECISION ALPHA,BETA
             #      INTEGER INCX,INCY,LDA,M,N
             #      CHARACTER TRANS
             #*     .. Array Arguments ..
             #      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
        function gemv!(trans::Char, m::Int, n::Int, alpha::($elty), A::Ptr{$elty}, stA::Int,
                        X::Ptr{$elty}, incX::Int, beta::($elty), Y::Ptr{$elty}, incY::Int)
            ccall((@blasfunc($fname), libblas), Void,
                (Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                 Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
                 Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}),
                 &trans, &m, &n, &alpha,
                 A, &stA, X, &incX,
                 &beta, Y, &incY)
            Y
        end

        function gemv!(trans::Char, alpha::($elty), A, X, beta::($elty), Y)
            m,n = size(A,1),size(A,2)
            if trans == 'N' && (length(X) != n || length(Y) != m)
                throw(DimensionMismatch("A has dimensions $(size(A)), X has length $(length(X)) and Y has length $(length(Y))"))
            elseif trans == 'C' && (length(X) != m || length(Y) != n)
                throw(DimensionMismatch("A' has dimensions $n, $m, X has length $(length(X)) and Y has length $(length(Y))"))
            elseif trans == 'T' && (length(X) != m || length(Y) != n)
                throw(DimensionMismatch("A.' has dimensions $n, $m, X has length $(length(X)) and Y has length $(length(Y))"))
            end
            gemv!(trans, m, n, alpha, pointer(A), max(1,stride(A,2)),
                    pointer(X), stride(X,1), beta, pointer(Y), stride(Y,1))
            Y
        end
    end
end

## Level 3


# We add some missing Ptr variants

gemm!(x...) = BLAS.gemm!(x...)

for (gemm, elty) in
        ((:dgemm_,:Float64),
         (:sgemm_,:Float32),
         (:zgemm_,:Complex128),
         (:cgemm_,:Complex64))
    @eval begin
             # SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
             # *     .. Scalar Arguments ..
             #       DOUBLE PRECISION ALPHA,BETA
             #       INTEGER K,LDA,LDB,LDC,M,N
             #       CHARACTER TRANSA,TRANSB
             # *     .. Array Arguments ..
             #       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
        function gemm!(transA::Char, transB::Char, m::Int, n::Int, k::Int, alpha::($elty),
                        A::Ptr{$elty}, stA::Int, B::Ptr{$elty}, stB::Int,
                        beta::($elty), C::Ptr{$elty},stC::Int)
            ccall((@blasfunc($gemm), libblas), Void,
                (Ptr{UInt8}, Ptr{UInt8}, Ptr{BlasInt}, Ptr{BlasInt},
                 Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt},
                 Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty},
                 Ptr{BlasInt}),
                 &transA, &transB, &m, &n,
                 &k, &alpha, A, &stA,
                 B, &stB, &beta, C,
                 &stC)
            C
        end

        function gemm!(transA::Char, transB::Char, alpha::($elty), A, B, beta::($elty), C)
#           if any([stride(A,1), stride(B,1), stride(C,1)] .!= 1)
#               error("gemm!: BLAS module requires contiguous matrix columns")
#           end  # should this be checked on every call?
            m = size(A, transA == 'N' ? 1 : 2)
            ka = size(A, transA == 'N' ? 2 : 1)
            kb = size(B, transB == 'N' ? 1 : 2)
            n = size(B, transB == 'N' ? 2 : 1)
            if ka != kb || m != size(C,1) || n != size(C,2)
                throw(DimensionMismatch("A has size ($m,$ka), B has size ($kb,$n), C has size $(size(C))"))
            end
            gemm!(transA, transB, m, n, ka, alpha,
                    pointer(A), max(1,stride(A,2)),
                    pointer(B), max(1,stride(B,2)), beta,
                    pointer(C), max(1,stride(C,2)))
            C
        end
    end
end
