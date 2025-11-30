module ApproxFun_Runtests

using ApproxFun
using ParallelTestRunner

# Start with autodiscovered tests
testsuite = find_tests(pwd())

# Parse arguments
args = parse_args(ARGS)

if filter_tests!(testsuite, args)
    # Remove tests that shouldn't run on Windows
    delete!(testsuite, "testutils")
end

const init_code = quote
    using ApproxFun
    using Random
    using Test
    using LinearAlgebra
    using SpecialFunctions
    using ApproxFunBase
    using ApproxFunBase: blocklengths, ∞, blockbandwidths, subblockbandwidths
    using ApproxFunBase.TestUtils: testbandedblockbandedoperator

    include(joinpath(@__DIR__, "testutils.jl"))
end

runtests(ApproxFun, args; testsuite, init_code)

end # module
