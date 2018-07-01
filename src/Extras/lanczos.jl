export lanczos


# this finds the OPs and recurrence for a
function lanczos(w,N)
    x = Fun(identity,space(w))

    f1=Fun(1/sqrt(sum(w)),space(x))

    P = Array{Fun}(undef, N + 1)
    β = Array{eltype(w)}(undef, N)
    γ = Array{eltype(w)}(undef, N)

    P[1] = f1

    v = x*P[1]
    β[1] = sum(w*v*P[1])

    v = v - β[1]*P[1]
    γ[1] = sqrt(sum(w*v^2))

    P[2] = v/γ[1]

    for k = 2:N
        v = x*P[k] - γ[k-1]*P[k-1]
        β[k] = sum(w*v*P[k])
        v = v - β[k]*P[k]
        γ[k] = sqrt(sum(w*v^2))
        P[k+1] = v/γ[k]
    end

    P,β,γ
end
