# Bootstraps shared by the estimators: did's multiplier bootstrap and the sup-t critical value

const TINY = sqrt(eps()) * 10

"""R's quantile of type 1: the smallest value with at least a share `q` of `v` at or below it."""
type1_quantile(v, q) = sort(v)[ceil(Int, q * length(v))]

"""The interquartile range of the draws `b` over the normal's, did's bootstrap standard deviation."""
normal_spread(b) = (type1_quantile(b, 0.75) - type1_quantile(b, 0.25)) / (quantile(Normal(), 0.75) - quantile(Normal(), 0.25))

"""Whether the draws `b` are missing or all but zero, as did's."""
degenerate(b) = any(isnan, b) || sum(abs2, b) <= TINY

"""
did's multiplier bootstrap (mboot): the clusters' sums of the influence functions (`S`, a row a cluster and a column
an estimate) each weighed by +1 or -1, a row a draw, scaled as did's.
"""
multiplier_draws(rng, S, draws) = (rand(rng, (-1.0, 1.0), draws, size(S, 1)) * S) ./ sqrt(size(S, 1))

"""
did's simultaneous critical value from the draws `B` of `multiplier_draws`: the 95% quantile of the largest |t| over
the columns that are not degenerate, each over its spread; nothing when there are none.
"""
function multiplier_critical(B)
    B = B[:, [k for k in axes(B, 2) if !degenerate(B[:, k])]]
    size(B, 2) == 0 && return nothing
    sigma = [(x = normal_spread(B[:, k]); x <= TINY ? NaN : x) for k in axes(B, 2)]
    largest = [maximum((x = abs(B[d, k] / sigma[k]); isnan(x) ? -Inf : x) for k in axes(B, 2)) for d in axes(B, 1)]
    filter!(isfinite, largest)
    isempty(largest) ? nothing : type1_quantile(largest, 0.95)
end

"""
The sup-t critical value of Montiel Olea and Plagborg-Moller (2019) for estimates with covariance `V`: the 95% quantile
of the largest |t| over those with a variance, in `draws` draws from their correlations, normal, or Student's t with
`df` degrees of freedom (a normal draw over the square root of a chi-squared draw over df) unless `df` is nothing;
nothing with fewer than two.
"""
function supt_critical(rng, V, draws, df)
    keep = [j for j in axes(V, 1) if !isnan(V[j, j]) && V[j, j] > 0]
    length(keep) <= 1 && return nothing
    sd = sqrt.([V[j, j] for j in keep])
    roots = eigen(Symmetric(V[keep, keep] ./ (sd * sd')))
    Z = randn(rng, draws, length(keep)) * (roots.vectors * Diagonal(sqrt.(max.(roots.values, 0.0))))'
    df === nothing || (Z ./= sqrt.(rand(rng, Chisq(df), draws) ./ df))
    quantile(vec(maximum(abs, Z; dims = 2)), 0.95)
end

"""
    wild_test(R, y, X, units, clusters, w; draws, rng, auxiliary, interval)

The wild cluster bootstrap test (Roodman, MacKinnon, Nielsen and Webb 2019) that `R` times the coefficients of the
weighted regression of `y` on `X` and unit effects is zero, the null imposed, the clusters' residuals each weighed by
`auxiliary` weights ("rademacher", every pattern of signs when there are no more than `draws`, or "webb"): its p-value,
with one row the 95% interval by inverting the test, and the draws. Defined when WildBootTests.jl is loaded.
"""
function wild_test end
