using RandomCorrelationMatrices
using Base.Test

srand(12345)
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