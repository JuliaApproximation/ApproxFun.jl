function deflationcheck(M::Matrix, Λ::Vector, tol::Real, verbose::Bool)
    b, k = size(M)
    @assert k == length(Λ)
    maxΛ = maximum(abs, Λ)
    minΛ = minimum(abs, Λ)
    j = 1
    maxΛ2j = (abs(Λ[j]) + minΛ)*(abs(Λ[j]) + maxΛ)
    FroM2 = colnorm2(M, j)
    while FroM2 ≤ tol*maxΛ2j && j < k
        j += 1
        maxΛ2j = max(maxΛ2j, (abs(Λ[j]) + minΛ)*(abs(Λ[j]) + maxΛ))
        FroM2 += colnorm2(M, j)
        verbose && println("    j = ", j, ", ||M_{1:b,1:j}||_F^2 = ", FroM2, ", ϵ (|λ_j| + |λ_min|)(|λ_j| + |λ_max|)  = ", tol*maxΛ2j, ".")
    end
    j-1
end

function deflationcheck(M::Matrix, N::Matrix, Λ::Vector, tol::Real, verbose::Bool)
    b, k = size(M)
    β, κ = size(N)
    @assert b == β && k == κ
    @assert k == length(Λ)
    maxΛ = maximum(abs, Λ)
    minΛ = minimum(abs, Λ)
    j = 1
    maxΛj = abs2(Λ[j])
    maxΛ2j = (abs(Λ[j]) + minΛ)*(abs(Λ[j]) + maxΛ)
    FroMN2 = colnorm2(M, j) + maxΛj*colnorm2(N, j)
    while FroMN2 ≤ tol*maxΛ2j && j < k
        j += 1
        maxΛj = max(maxΛj, abs2(Λ[j]))
        maxΛ2j = max(maxΛ2j, (abs(Λ[j]) + minΛ)*(abs(Λ[j]) + maxΛ))
        FroMN2 += colnorm2(M, j) + maxΛj*colnorm2(N, j)
        verbose && println("    j = ", j, ", ||M_{1:b,1:j}||_F^2 + |λ_j|^2||N_{1:b,1:j}||_F^2 = ", FroMN2, ", ϵ (|λ_j| + |λ_min|)(|λ_j| + |λ_max|)  = ", tol*maxΛ2j, ".")
    end
    j-1
end

function colnorm2(M::Matrix, J::Int)
    ret = zero(real(eltype(M)))
    @inbounds for i = 1:size(M, 1)
        ret += abs2(M[i, J])
    end
    ret
end
