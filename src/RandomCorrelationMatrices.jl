module RandomCorrelationMatrices

    using Distributions
    
    # Based on
    # Lewandowski, Kurowicka, Joe
    # "Generating random correlation matrices based on vines
    #   and extended onion method"
    # Journal of Multivariate Analysis 100 (2009)
    # doi:10.1016/j.jmva.2009.04.008
    export randcormatrix
    function randcormatrix(d, η)
        β = η + (d-2)/2
        u = rand(Beta(β,β))
        r₁₂ = 2*u - 1
        r = [   1 r₁₂ ;
              r₁₂   1 ]
        # In published paper this index is n = 2:d-1
        # n is never mentioned again, but a mysterious k is
        # We'll thus use k
        # Similar approach was taken here:
        # http://stats.stackexchange.com/a/125017/58921
        for k in 2:d-1
            β = β - 1/2
            y = rand(Beta(k/2,β))
            u = randn(k)
            u /= norm(u)
            w = sqrt(y)*u
            A = chol(r, :L)
            z = A*w
            r = [ r  z ;
                  z' 1 ]
        end
        return r
    end


    # Use randcormatrix + a desired vector of standard deviations
    # for each term to generate a random covariance matrix
    export randcovmatrix
    function randcovmatrix(d, η, σ)
        length(σ) != d && throw(DimensionMismatch("length(sigma) doesn't match d"))
        r = randcormatrix(d, η)
        Σ = zeros(d,d)
        @inbounds for i in 1:d, j in 1:d
            Σ[i,j] = r[i,j] * σ[i] * σ[j]
        end
        return Σ
    end

end # module
