using Documenter, ApproxFun

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/ApproxFun/ApproxFun.jl.git",
    latest = "development",
    julia  = "0.5",
    osname = "linux")
