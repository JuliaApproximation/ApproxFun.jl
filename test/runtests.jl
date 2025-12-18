module ApproxFun_Runtests

using ApproxFun
using ParallelTestRunner

# Start with autodiscovered tests
testsuite = find_tests(pwd())

if "--downstream_integration_test" in ARGS
    delete!(testsuite, "test_aqua")
end
filtered_args = filter(!=("--downstream_integration_test"), ARGS)

# Parse arguments
args = parse_args(filtered_args)

if filter_tests!(testsuite, args)
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
