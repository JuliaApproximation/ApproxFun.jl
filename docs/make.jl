using Documenter, ApproxFun

makedocs(modules=[ApproxFun],
			doctest = true,
			clean = true,
			format = :html,
			sitename = "ApproxFun.jl",
			authors = "Sheehan Olver",
			pages = Any[
					"Home" => "index.md",
					"Usage" => Any[
						"Constructors" => "usage/constructors.md"
					],
					"Frequently Asked Questions" => "faq.md",
					"Library" => "library.md"
				]
			)


deploydocs(
    repo   = "github.com/ApproxFun/ApproxFun.jl.git",
    latest = "development",
    julia  = "0.5",
    osname = "linux",
    target = "build",
    deps   = nothing,
    make   = nothing
    )
