using Documenter
DocMeta.setdocmeta!(ApproxFunBase, :DocTestSetup, :(using ApproxFun); recursive=true)
DocMeta.setdocmeta!(ApproxFun, :DocTestSetup, :(using ApproxFun); recursive=true)

@testset "doctests" begin
    doctest(ApproxFun)
    doctest(ApproxFunBase, manual=false)
end