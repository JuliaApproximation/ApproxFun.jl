

Base.eigvals(A::BandedOperator,n::Integer)=eigvals(full(A[1:n,1:n]),full(Conversion(domainspace(A),rangespace(A))[1:n,1:n]))