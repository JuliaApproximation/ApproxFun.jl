## Testing
# These routines are for the unit tests

using Base.Test


## Supports @test_approx_eq


Base.Test.approx_full(f::Fun) = f



## Spaces Tests


function testtransforms(S::Space;minpoints=1,invertibletransform=true)
    # transform tests
    v = rand(max(minpoints,min(100,ApproxFun.dimension(S))))
    plan = plan_transform(S,v)
    @test transform(S,v)  == plan*v

    iplan = plan_itransform(S,v)
    @test itransform(S,v)  == iplan*v

    if invertibletransform
        for k=max(1,minpoints):min(5,dimension(S))
            v = [zeros(k-1);1.0]
            @test transform(S,itransform(S,v)) ≈ v
        end

        @test transform(S,itransform(S,v)) ≈ v
        @test itransform(S,transform(S,v)) ≈ v
    end
end

function testcalculus(S::Space;haslineintegral=true,hasintegral=true)
    for k=1:min(5,dimension(S))
        v = [zeros(k-1);1.0]
        f = Fun(S,v)
        @test abs(DefiniteIntegral()*f-sum(f)) < 100eps()
        if haslineintegral
            @test DefiniteLineIntegral()*f ≈ linesum(f)
        end
        @test norm(Derivative()*f-f') < 100eps()
        if hasintegral
            @test norm(differentiate(integrate(f))-f) < 100eps()
            @test norm(differentiate(cumsum(f))-f) < 100eps()
            @test norm(first(cumsum(f))) < 100eps()
        end
    end
end

function testmultiplication(spa,spb)
    for k=1:10
        a = Fun(spa,[zeros(k-1);1.])
        M = Multiplication(a,spb)
        pts = ApproxFun.checkpoints(rangespace(M))
        for j=1:10
            b = Fun(spb,[zeros(j-1);1.])
            @test (M*b).(pts) ≈ a.(pts).*b.(pts)
        end
    end
end

function testspace(S::Space;
                    minpoints=1,invertibletransform=true,haslineintegral=true,hasintegral=true,
                    dualspace=S)
    testtransforms(S;minpoints=minpoints,invertibletransform=invertibletransform)
    testcalculus(S;haslineintegral=haslineintegral,hasintegral=hasintegral)
    if dualspace ≠ nothing
        testmultiplication(dualspace,S)
    end
end





## Operator Tests

function backend_testfunctional(A)
    @test rowstart(A,1) ≥ 1
    @test colstop(A,1) ≤ 1
    @test bandwidth(A,1) ≤ 0
    @test blockbandwidth(A,1) ≤ 0

    B=A[1:10]
    @test eltype(B) == eltype(A)
    for k=1:5
        @test B[k] ≈ A[k]
        @test isa(A[k],eltype(A))
    end
    @test B ≈ A[1,1:10]
    @test B.' ≈ A[1:1,1:10]
    @test B[3:10] ≈ A[3:10]
    @test B ≈ [A[k] for k=1:10]



    co=cache(A)
    @test co[1:10] ≈ A[1:10]
    @test co[1:10] ≈ A[1:10]
    @test co[20:30] ≈ A[1:30][20:30] ≈ A[20:30]
end

# Check that the tests pass after conversion as well
function testfunctional{T<:Real}(A::Operator{T})
    backend_testfunctional(A)
    backend_testfunctional(Operator{Float64}(A))
    backend_testfunctional(Operator{Float32}(A))
    backend_testfunctional(Operator{Complex128}(A))
end

function testfunctional{T<:Complex}(A::Operator{T})
    backend_testfunctional(A)
    backend_testfunctional(Operator{Complex64}(A))
    backend_testfunctional(Operator{Complex128}(A))
end

function backend_testinfoperator(A)
    @test isinf(size(A,1))
    @test isinf(size(A,2))
    B=A[1:5,1:5]
    eltype(B) == eltype(A)

    for k=1:5,j=1:5
        @test B[k,j] ≈ A[k,j]
        @test isa(A[k,j],eltype(A))
    end

    @test A[1:5,1:5][2:5,1:5] ≈ A[2:5,1:5]
    @test A[1:5,2:5] ≈ A[1:5,1:5][:,2:end]
    @test A[1:10,1:10][5:10,5:10] ≈ [A[k,j] for k=5:10,j=5:10]
    @test A[1:10,1:10][5:10,5:10] ≈ A[5:10,5:10]
    @test A[1:30,1:30][20:30,20:30] ≈ A[20:30,20:30]

    @test A[Block(1):Block(3),Block(1):Block(3)] ≈ A[blockstart(rangespace(A),1):blockstop(rangespace(A),3),blockstart(domainspace(A),1):blockstop(domainspace(A),3)]
    @test A[Block(3):Block(4),Block(2):Block(4)] ≈ A[blockstart(rangespace(A),3):blockstop(rangespace(A),4),blockstart(domainspace(A),2):blockstop(domainspace(A),4)]

    for k=1:10
        @test isfinite(colstart(A,k)) && colstart(A,k) > 0
        @test isfinite(rowstart(A,k)) && colstart(A,k) > 0
    end

    co=cache(A)
    @test co[1:10,1:10] ≈ A[1:10,1:10]
    @test co[1:10,1:10] ≈ A[1:10,1:10]
    @test co[20:30,20:30] ≈ A[1:30,1:30][20:30,20:30]

    let C=cache(A)
        resizedata!(C,5,35)
        resizedata!(C,10,35)
        @test C.data[1:10,1:C.datasize[2]] ≈ A[1:10,1:C.datasize[2]]
    end
end

# Check that the tests pass after conversion as well
function testinfoperator{T<:Real}(A::Operator{T})
    backend_testinfoperator(A)
    backend_testinfoperator(Operator{Float64}(A))
    backend_testinfoperator(Operator{Float32}(A))
    backend_testinfoperator(Operator{Complex128}(A))
end

function testinfoperator{T<:Complex}(A::Operator{T})
    backend_testinfoperator(A)
    backend_testinfoperator(Operator{Complex64}(A))
    backend_testinfoperator(Operator{Complex128}(A))
end

function testraggedbelowoperator(A)
    @test israggedbelow(A)
    for k=1:20
        @test isfinite(colstop(A,k))
    end
    testinfoperator(A)
end

function testbandedbelowoperator(A)
    @test isbandedbelow(A)
    @test isfinite(bandwidth(A,1))
    testraggedbelowoperator(A)

    for k=1:10
        @test colstop(A,k) ≤ k + bandwidth(A,1)
    end
end


function testalmostbandedoperator(A)
    testbandedbelowoperator(A)
end

function testbandedoperator(A)
    @test isbanded(A)
    @test isfinite(bandwidth(A,2))
    testalmostbandedoperator(A)
    for k=1:10
        @test rowstop(A,k) ≤ k + bandwidth(A,2)
    end

    @test isa(A[1:10,1:10],BandedMatrix)
end


function testblockbandedoperator(A)
    @test isblockbanded(A)
    testraggedbelowoperator(A)
    @test isfinite(blockbandwidth(A,2))
    @test isfinite(blockbandwidth(A,1))

    for K=1:10
        @test K - blockbandwidth(A,2) ≤ blockcolstop(A,K).K ≤ K + blockbandwidth(A,1) < ∞
        @test K - blockbandwidth(A,1) ≤ blockrowstop(A,K).K ≤ K + blockbandwidth(A,2) < ∞
    end
end

function testbandedblockbandedoperator(A)
    @test isbandedblockbanded(A)
    testblockbandedoperator(A)
    @test isfinite(subblockbandwidth(A,1))
    @test isfinite(subblockbandwidth(A,2))

    @test isa(A[1:10,1:10],BandedBlockBandedMatrix)
end
