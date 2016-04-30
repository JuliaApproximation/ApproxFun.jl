


pad!(A::BandedMatrix,n,::Colon)=pad!(A,n,n+A.u)  # Default is to get all columns
