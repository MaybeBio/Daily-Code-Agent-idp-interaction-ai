using Droplet

function euler_maruyama(x0, kappa, gamma, D, dt, n; seed=nothing)
    x = Vector{Float64}(undef, n)
    x[1] = x0
    xi = randn(n-1)
    for i in 1:n-1
        x[i+1] = x[i] - (kappa/gamma)*x[i]*dt + sqrt(2*D*dt)*xi[i]
    end
    return x
end

let
    start_time = time()

    x0 = 0.0
    kappa = 1.0
    gamma = 1.0
    D = 1.0
    dt = 0.01 * (gamma / kappa)
    n = 10_000_000

    euler_maruyama(x0, kappa, gamma, D, dt, n)

    end_time = time()

    print(end_time - start_time)

end