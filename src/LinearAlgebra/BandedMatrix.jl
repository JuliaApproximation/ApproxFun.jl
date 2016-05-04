


pad!(A::BandedMatrix,n,::Colon) = pad!(A,n,n+A.u)  # Default is to get all columns


columnrange(A,row::Integer) = max(1,row+bandinds(A,1)):row+bandinds(A,2)
