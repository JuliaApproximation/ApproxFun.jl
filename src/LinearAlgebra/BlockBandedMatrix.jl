immutable BlockBandedMatrix{T,RI,CI} <: AbstractMatrix{T}
    data::Vector{T}  # the entries
    rows::RI
    columns::CI

    l::Int # lower bandwidth ≥0
    u::Int # upper bandwidth ≥0

    columnlookup::Vector{Int}  # translates column first data entry
    firstrow::Vector{Int}     # gives first row index in

    function BlockBandedMatrix(dat::Vector{T},rws::RI,cls::CI,l::Int,u::Int)
        n = 0
        for J=1:length(columns),K=max(1,J-u):min(J+l,length(rows))
            n += rows[K]*columns[J]
        end

        if n != length(dat)
            error("Data must be consistant size with number of non-zeros")
        end

        m = 1
        curcol = 1

        columnlookup = Array(Int,0)
        firstrow = Array(Int,0)
        for J=1:length(columns)
            sm = sum(rows[max(1,J-u):min(J+l,length(rows))])
            nrws = sum(rows[max(1,J-u):J-1])
            for j=1:columns[J]
                push!(columnlookup,m)
                push!(firstrow,curcol-j-nrws+1)

                m += sm
                curcol += 1
            end
        end

        new(dat,rws,cls,l,u,columnlookup,firstrow)
    end
end

BlockBandedMatrix(dat::Vector,rws,cls,l,u) =
    BlockBandedMatrix{eltype(dat),typeof(rws),typeof(cls)}(dat,rws,cls,l,u)

Base.size(B::BlockBandedMatrix) = sum(B.rows),sum(B.columns)
Base.linearindexing{BBM<:BlockBandedMatrix}(::Type{BBM}) = Base.LinearSlow()
function Base.getindex(B::BlockBandedMatrix,k::Integer,j::Integer)
    if k > size(B,1) || j > size(B,2)
        throw(BoundsError())
    end

    fr = B.firstrow[j]
    nrws = j != size(B,2) ?
            B.columnlookup[j+1] - B.columnlookup[j] :
            length(B.data) - B.columnlookup[j] + 1

    if k < fr || k > fr + nrws - 1
        zero(eltype(B))
    else
        B.data[B.columnlookup[j] + k - fr]
    end
end

function Base.setindex!(B::BlockBandedMatrix,v,k::Integer,j::Integer)
    if k > size(B,1) || j > size(B,2)
        throw(BoundsError())
    end

    fr = B.firstrow[j]
    nrws = j != size(B,2) ?
            B.columnlookup[j+1] - B.columnlookup[j] :
            length(B.data) - B.columnlookup[j] + 1

    if k < fr || k > fr + nrws - 1
        if v!=0
            error("Cannot change entries outside block bands")
        end
    else
        B.data[B.columnlookup[j] + k - fr] = v
    end
end


for (op,bop) in ((:(Base.rand),:bbrand),(:(Base.zeros),:bbzeros),(:(Base.ones),:bbones))
    @eval begin
        function $bop{T}(::Type{T},rows,columns,l,u)
            n = 0
            for J=1:length(columns),K=max(1,J-u):min((J+l),length(rows))
                n += rows[K]*columns[J]
            end
            BlockBandedMatrix($op(T,n),rows,columns,l,u)
        end
        $bop(rows,columns,l,u) = $bop(Float64,rows,columns,l,u)
    end
end




rows=1:50
columns=1:50
l=u=1

B=bbrand(rows,columns,l,u)


sz=sizeof(eltype(B))

v=pointer(B.data)
α=1.0
Z=bbzeros(rows,columns,l,u+l)
    w=pointer(Z.data)
    for k=1:4
        N=B.columnlookup[k+1]-B.columnlookup[k]
        BLAS.axpy!(N,α,v+(B.columnlookup[k]-1)*sz,1,w+(B.columnlookup[k]-1)*sz,1)
    end


B
unsafe_load(v+1*sz)
Z
using SO
B

B[10,3]=5

B

v=B.data[1:3]
norm(v)






B.columnlookup[2]





B.columnlookup

B.firstrow

size(B)

k=j=1
    B.data[B.columnlookup[j] + k-1]



B.data



columnlookup = Array(Int,0)
    m=1
    push!(columnlookup,m)  # the first column begins at index 1
    for J=1:length(columns)
        sm = sum(rows[max(1,J-u):min(J+l,length(rows))])
        for j=1:columns[J]
            m += sm
            push!(columnlookup,m)
        end
    end



columns[2]

columnlookup



n = 0
    for J=1:length(columns),K=max(1,J-u):(J+l)
        n += rows[K]*columns[J]
    end


    B=BlockBandedMatrix(rand(n),1,1,rows,columns)




B.data

size(B)



T=eltype(data)
@time pointer_to_array(pointer(data)+sizeof(T)*2,(2,2))
