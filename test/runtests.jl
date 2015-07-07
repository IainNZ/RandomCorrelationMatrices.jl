using RandomCorrelationMatrices
using Base.Test

srand(12345)

println("Random correlation matrix")
n = 4
for η in [2,8,32,128]
    r = randcormatrix(n,η)
    @show η
    for i in 1:n
        for j in 1:n
            @printf("%5.2f ", r[i,j])
        end
        println()
    end

    ranges = Float64[]
    for s in 1:1000
        r = randcormatrix(n,η)
        r0 = r - diagm(diag(r))
        push!(ranges, maximum(vec(r0)) - minimum(vec(r0)))
    end
    @show mean(ranges)
    @show median(ranges)
end


println("\nRandom covariance matrix")
σ = [100, 200, 300, 400, 500]
@show σ
Σ = randcovmatrix(5, 10, σ)
L = chol(Σ,:L)
samples = 100_000
results = zeros(5,samples)
for i in 1:samples
    results[:,i] = L*randn(5)
end
@show round(mean(results,2),1)
@show round(std(results,2),1)