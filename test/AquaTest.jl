using Aqua
@testset "Project quality" begin
    Aqua.test_all(ApproxFun, ambiguities=false, undefined_exports=false,
        stale_deps=(; ignore=[:ApproxFunBaseTest]), piracies = false)
end
