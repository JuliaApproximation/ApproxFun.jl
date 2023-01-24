macro verbose(ex)
    head = ex.head
    args = ex.args
    @assert args[1] == Symbol("@testset")
    name = args[3] isa String ? args[3] : nothing
    if VERSION >= v"1.8"
        insert!(args, 3, Expr(:(=), :verbose, true))
    end
    Expr(head, args...)
end
