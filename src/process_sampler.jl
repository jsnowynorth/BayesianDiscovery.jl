

"""
hpd(x::Vector; p = 0.95)

Computes the highest posterior density of a vector x
"""
# function hpd(x::Vector; p=0.95)
#     p1 = round((1 - p) / 2, digits=5)
#     p2 = round(p + (1 - p) / 2, digits=5)

#     # dens = kde(x, boundary=extrema(x))
#     # # dens = kde(x)
#     # Δn = dens.x[2] - dens.x[1]
#     # inds = (cumsum(dens.density * Δn) .> p1) .& (cumsum(dens.density * Δn) .< p2)

#     # # lower = dens.x[inds][1]
#     # # upper = dens.x[inds][end]

#     # if sum(inds) < 2
#     #     lower = 0
#     #     upper = 0
#     # else
#     #     lower = dens.x[inds][1]
#     #     upper = dens.x[inds][end]
#     # end

#     df = DataFrame(x=sort(x), y=range(0, 1, length=length(x)))
#     lower = df[:, :x][(df[:, :y].>p1).&(df[:, :y].<p2)][1]
#     upper = df[:, :x][(df[:, :y].>p1).&(df[:, :y].<p2)][end]

#     return lower, upper
# end

# https://rpubs.com/tjmahr/hpdi
function hpd(x; p=0.95)
    if isempty(x)
        return [0,0]
    else
        x_sorted = sort(x)
        n = length(x)
        num_to_keep = Int(ceil(p * n))
        num_to_drop = n - num_to_keep

        possible_starts = range(1, num_to_drop + 1, step=1)
        # Just count down from the other end
        possible_ends = range(n - num_to_drop, n, step=1)

        # Find smallest interval
        span = x_sorted[possible_ends] - x_sorted[possible_starts]
        edge = argmin(span)
        edges = [possible_starts[edge], possible_ends[edge]]

        return x_sorted[edges]
    end
    
end



struct Posterior_Summary

    M
    M_hpd
    gamma
    ΣZ
    ΣU
    π

end

function posterior_summary(model::Model, pars::Pars, posterior::Posterior)

    M_mu = [round(mean(posterior.M[n, d, posterior.gamma[n,d,:] .!= 0]), digits = 3) for n in 1:model.N, d in 1:model.D]
    # M_mu = mean(posterior.M, dims=3)
    gamma_mu = mean(posterior.gamma, dims=3)
    ΣZ_mu = mean(posterior.ΣZ, dims=3)
    ΣU_mu = mean(posterior.ΣU, dims=3)
    π_mu = mean(posterior.π, dims=2)

    # M_hpd = hcat(collect.([hpd(posterior.M[n, d, :]) for n in 1:model.N, d in 1:model.D])...)
    M_hpd = hcat(collect.([hpd(posterior.M[n, d, posterior.gamma[n,d,:] .!= 0]) for n in 1:model.N, d in 1:model.D])...)


    post_sum = Posterior_Summary(M_mu, M_hpd, gamma_mu, ΣZ_mu, ΣU_mu, π_mu)

    return post_sum

end

function posterior_surface(model::Model, pars::Pars, posterior::Posterior)

    Ψ = model.Basis.SpaceDerivative[model.Basis.SpaceNames[1]]
    Φ = model.Basis.TimeDerivative[model.Basis.TimeNames[1]]

    nSpace = size(Ψ, 2)
    nTime = size(Φ, 2)
    nits = size(posterior.A, 3)

    if model.N == 1
        postA = cat([tensor_mult(fold3(posterior.A[:, :, i], nSpace, nTime, model.N), Ψ, Φ, model.Θ) for i in 1:nits]..., dims=3)
        return mean(postA, dims=3), std(postA, dims=3)
    else
        postA = cat([tensor_mult(fold3(posterior.A[:, :, i], nSpace, nTime, model.N), Ψ, Φ, model.Θ) for i in 1:nits]..., dims=4)
        return dropdims(mean(postA, dims=4), dims=4), dropdims(std(postA, dims=4), dims=4)
    end

    # postA = cat([tensor_mult(fold3(posterior.A[:, :, i], nSpace, nTime, model.N), Ψ, Φ, model.Θ) for i in 1:nits]..., dims=3)

    # return mean(postA, dims=3), std(postA, dims=3)

end


struct Equation

    mean
    lower
    upper
    percentile
    inclusion_probability

end

Base.show(io::IO, equation::Equation) =
    print(io, "Equation\n",
        " ├─── Mean:      ", equation.mean, '\n',
        " ├─── Lower HPD: ", equation.lower, '\n',
        " ├─── Upper HPD: ", equation.upper, '\n',
        " ├─── Percentile: ", equation.percentile, '\n',
        " └─── Inclusion Probability: ", equation.inclusion_probability)
#

function print_equation(sys_names, model::Model, pars::Pars, posterior::Posterior; cutoff_prob=0.5, p=0.95)

    # M_est = round.(mean(posterior.M, dims=3), digits=3)
    # M_hpd = hcat(collect.([hpd(posterior.M[n, d, :]; p=p) for n in 1:model.N, d in 1:model.D])...)
    M_est = [round(mean(posterior.M[n, d, posterior.gamma[n,d,:] .!= 0]), digits = 3) for n in 1:model.N, d in 1:model.D]
    M_hpd = hcat(collect.([hpd(posterior.M[n, d, posterior.gamma[n,d,:] .!= 0]) for n in 1:model.N, d in 1:model.D])...)
    gamma_mu = mean(posterior.gamma, dims=3)
    # a = size(M_est)[1]
    # b = size(M_est)[2]

    lower = reshape(round.(M_hpd[1, :], digits=3), model.N, model.D)
    upper = reshape(round.(M_hpd[2, :], digits=3), model.N, model.D)

    par_names = reshape(split(model.function_names, ", "), 1, model.D)

    # if model.X === nothing
    #     nms = replace((@code_string model.Λ(model.Z[:, 1], [ones(model.N) for i in 1:(length(model.Λnames)-1)])).string[findfirst("return", (@code_string model.Λ(model.Z[:, 1], [ones(model.N) for i in 1:(length(model.Λnames)-1)])).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "")
    #     par_names = reshape(split(nms, ", "), 1, model.D)
    # else
    #     nms = replace((@code_string model.Λ(model.Z[:, 1], [ones(model.N) for i in 1:(length(model.Λnames)-1)], model.X[:, 1])).string[findfirst("return", (@code_string model.Λ(model.Z[:, 1], [ones(model.N) for i in 1:(length(model.Λnames)-1)], model.X[:, 1])).string)[end]:end], "\n    " => " ", "\n" => "", "n" => "", "end" => "", "[" => "", "]" => "")
    #     par_names = vcat(reshape(split(nms, ", "), 1, model.D))
    # end


    included = gamma_mu .> cutoff_prob
    final_eqs = copy(sys_names)
    final_lower = copy(sys_names)
    final_upper = copy(sys_names)

    n_eqs = size(final_eqs, 1)

    for i in 1:n_eqs
        found = map((x, y) -> string(x, ' ', y), M_est[i, included[i, :]], par_names[included[i, :]])
        no_sign = replace.(found, "-" => "")
        signs = sign.(M_est[i, included[i, :]])
        signs = [(signs[i] == 1 ? "+" : "-") for i in 1:size(signs)[1]]
        if signs[1] == "+"
            replace(signs[1], "+" => "")
        end
        found = join(map((x, y) -> string(x, ' ', y), signs, no_sign), " ")

        final_eqs[i] = string(final_eqs[i], " = ", found)
    end

    for i in 1:n_eqs
        found = map((x, y) -> string(x, ' ', y), lower[i, included[i, :]], par_names[included[i, :]])
        no_sign = replace.(found, "-" => "")
        signs = sign.(lower[i, included[i, :]])
        signs = [(signs[i] == 1 ? "+" : "-") for i in 1:size(signs)[1]]
        if signs[1] == "+"
            replace(signs[1], "+" => "")
        end
        found = join(map((x, y) -> string(x, ' ', y), signs, no_sign), " ")

        final_lower[i] = string(final_lower[i], " = ", found)
    end

    for i in 1:n_eqs

        found = map((x, y) -> string(x, ' ', y), upper[i, included[i, :]], par_names[included[i, :]])
        no_sign = replace.(found, "-" => "")
        signs = sign.(upper[i, included[i, :]])
        signs = [(signs[i] == 1 ? "+" : "-") for i in 1:size(signs)[1]]
        if signs[1] == "+"
            replace(signs[1], "+" => "")
        end
        found = join(map((x, y) -> string(x, ' ', y), signs, no_sign), " ")

        final_upper[i] = string(final_upper[i], " = ", found)
    end

    equation = Equation(final_eqs, final_lower, final_upper, p, cutoff_prob)

    return equation
end
