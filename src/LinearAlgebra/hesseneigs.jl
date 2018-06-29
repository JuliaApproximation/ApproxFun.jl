


import LinearAlgebra.BLAS: BlasInt, @blasfunc


for (hseqr,elty) in ((:zhseqr_,:ComplexF64),)
    @eval function hesseneigvals(M::Matrix{$elty})
        if isempty(M)
            return $elty[]
        end
       A=vec(M)

        N=size(M,1)
        info  =0

        ilo = 1; ihi = N; ldh=N;ldz=N;lwork = N

        z=zero($elty)
        work  = Array{$elty}( N*N)
        w=Array{$elty}(N)

        Ec='E'
        Nc='N'
        ccall((@blasfunc($hseqr),LAPACK.liblapack),
            Nothing,
            (Ref{UInt8}, Ref{UInt8},
        Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, #A
        Ref{BlasInt}, Ptr{$elty}, Ref{$elty}, #z
        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}),
        Ec, Nc,
        N , ilo, ihi, A,
        ldh, w, z,
        ldz, work, lwork, info)
        w
    end
end

for (hseqr,elty) in ((:dhseqr_,:Float64),)
    @eval function hesseneigvals(M::Matrix{$elty})
        if isempty(M)
            return $elty[]
        end

        A=vec(M)

        N=size(M,1)
        info  =0

        ilo = 1; ihi = N; ldh=N;ldz=N;

        lwork = -1

        z=zero($elty)
        work  = Array{$elty}( 1)
        wr=Array{$elty}(N)
        wi=Array{$elty}(N)

        Ec='E'
        Nc='N'
        for i=1:2
            ccall((@blasfunc($hseqr),LAPACK.liblapack),
                Nothing,
                (Ref{UInt8},Ref{UInt8},
            Ref{BlasInt},Ref{BlasInt},Ref{BlasInt},Ptr{$elty}, #A
            Ref{BlasInt},Ptr{$elty},Ptr{$elty},Ref{$elty}, #z
            Ref{BlasInt},Ptr{$elty},Ref{BlasInt},Ref{BlasInt}),
            Ec, Nc,
            N , ilo, ihi, A,
            ldh, wr,wi, z,
            ldz, work, lwork, info)

            if lwork < 0
                lwork=Int(real(work[1]))
                work=Array{$elty}(lwork)
            end
        end

        wr+im*wi
    end
end
