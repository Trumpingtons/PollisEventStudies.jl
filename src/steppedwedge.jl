# The mixed model of a stepped-wedge cluster randomized trial (Hussey and Hughes 2007), with the effect by exposure time
# of Kenny et al. (2022): a port of the mixed-model part of R's steppedwedge 1.0.0 (Kenny, Arthur and Ouyang), fitted by
# MixedModels.jl

"""
    mixed_fit(data, fixed, groups, slope; family, weights)

The mixed model of `data`'s `y` on the `fixed` columns, with no intercept, a random intercept for each column of
`groups`, and with `slope` a random coefficient of it by `exposure`, fitted by REML when `family` is "gaussian" and by
the Laplace approximation otherwise: the fixed effects with their names and covariance, the standard deviations of the
random effects by grouping, that of the residual (nothing but for "gaussian"), and the slope's conditional modes and
variances by level of `exposure`. Defined when MixedModels.jl is loaded.
"""
function mixed_fit end

"""A stepped-wedge trial's result: the effect, its standard error and 95% interval, the effects by exposure time, the
standard deviations of the random effects and the intracluster correlation, the sample and the fitted model."""
struct SteppedWedgeResult
    # "IT", "ETI" or "TEH"
    effect::String
    # "TATE" or "PTE", at the exposure times `exposure` (a range, or a whole number)
    estimand::String
    exposure::Union{UnitRange{Int}, Int}
    estimate::Float64
    se::Float64
    interval::Tuple{Float64, Float64}
    # An exposure time a row: the effect at it and its standard error
    by_exposure::NamedTuple
    # The standard deviations of the random intercepts of the clusters, of the clusters at each time, of the people
    # (nothing without them) and of the effects by exposure time (nothing but under "TEH"), and of the residual (nothing
    # but for a gaussian outcome)
    sd::NamedTuple
    icc::Union{Float64, Nothing}
    # Observations, clusters, times and sequences (the distinct times clusters start the treatment); the family
    sample::NamedTuple
    # The MixedModels.jl model
    model::Any
end

"""
    stepped_wedge(table; cluster, time, treatment, outcome, trials = nothing, individual = nothing, effect = "IT",
        estimand = "TATE", exposure = nothing, family = "gaussian")

The effect of a treatment rolled out to clusters at random times, a stepped wedge: after `using MixedModels`, the mixed
model of Hussey and Hughes (2007), with a fixed effect for each time, a random intercept for each cluster and one for
each cluster at each time (so that the correlation of outcomes in a cluster is stronger at the same time than at
different times), as R's steppedwedge `analyze(method = "mixed")`.

`table` is any Tables.jl table, a row an observation (a person, or a cluster at a time); `cluster`, `time`, `treatment`
and `outcome` name its columns. `treatment` is 0 or 1, the same for all the observations of a cluster at a time. The
exposure time of a treated observation is its time less that of its cluster's first treatment, plus 1. With `trials`,
`outcome` counts the successes out of `trials`, and `family` is "binomial". With `family = "binomial"` and no `trials`,
`outcome` is 0 or 1. `individual` adds a random intercept for each person, seen at several times.

The `effect` of the treatment:

- "IT", immediate treatment: the same effect at every exposure time. When the effect grows (or fades) with exposure
  time, its estimate can be far from the effect averaged over the exposure times, even of the wrong sign (Kenny et al.
  2022);
- "ETI", exposure time indicators: an effect at each exposure time;
- "TEH", treatment effect heterogeneity: a random effect at each exposure time about a mean effect.

`estimand` is "TATE", the average of the effects over exposure times `exposure`, a range (all of them by default), or
"PTE", the effect at exposure time `exposure`, a whole number. Under "IT" both are the one effect; under "TEH" the TATE
is the mean effect, and a PTE the mean plus the exposure time's conditional mode, with the variance of both. The 95%
interval is the estimate plus or minus 1.96 standard errors, as steppedwedge's.

Gaussian models are fitted by REML, as R's lme4 `lmer`; binomial ones by the Laplace approximation, as `glmer`.
steppedwedge 1.0.0 takes the standard error of a PTE under "TEH" from the conditional variance of the exposure time
before the one asked for; here it is that of the exposure time asked for.

The result, a `SteppedWedgeResult`, has the estimate, its standard error and 95% interval, the effects by exposure time
(under "ETI" the fixed effects, under "TEH" the mean plus each conditional mode), the standard deviations of the random
effects and of the residual, the intracluster correlation (that of two people of the same cluster at the same time) of
a gaussian model, the sample and the fitted model.
"""
function stepped_wedge(table; cluster, time, treatment, outcome, trials = nothing, individual = nothing,
        effect::AbstractString = "IT", estimand::AbstractString = "TATE", exposure = nothing,
        family::AbstractString = "gaussian")
    Base.get_extension(@__MODULE__, :PollisEventStudiesMixedModelsExt) === nothing &&
        throw(ArgumentError("the stepped-wedge mixed model needs MixedModels.jl: using MixedModels"))
    effect in ("IT", "ETI", "TEH") || throw(ArgumentError("effect is \"IT\", \"ETI\" or \"TEH\", not \"$effect\""))
    estimand in ("TATE", "PTE") || throw(ArgumentError("estimand is \"TATE\" or \"PTE\", not \"$estimand\""))
    family in ("gaussian", "binomial") || throw(ArgumentError("family is \"gaussian\" or \"binomial\", not \"$family\""))
    trials === nothing || family == "binomial" ||
        throw(ArgumentError("an outcome of successes out of trials is binomial: take family = \"binomial\""))
    estimand == "PTE" && !(exposure isa Integer) &&
        throw(ArgumentError("the PTE is the effect at an exposure time: give it as exposure, a whole number"))
    estimand == "TATE" && !(exposure === nothing || exposure isa AbstractUnitRange{<:Integer}) &&
        throw(ArgumentError("the TATE averages the effects over a range of exposure times, such as exposure = 1:4"))

    columns = Tables.columns(table)
    column(name) = Tables.getcolumn(columns, Symbol(name))
    named = [cluster, time, treatment, outcome]
    trials === nothing || push!(named, trials)
    individual === nothing || push!(named, individual)
    needed = [column(name) for name in named]
    rows = [k for k in 1:length(needed[1]) if all(!ismissing(c[k]) for c in needed)]
    isempty(rows) && throw(ArgumentError("no observation has all of $(join(named, ", "))"))

    clusters = [string(column(cluster)[k]) for k in rows]
    times = [whole(column(time)[k]) for k in rows]
    treated = map(rows) do k
        d = column(treatment)[k]
        d in (0, 1) || throw(ArgumentError("treatment is 0 or 1, not $d"))
        Int(d)
    end
    cells = Dict{Tuple{String, Int}, Int}()
    for (c, t, d) in zip(clusters, times, treated)
        get!(cells, (c, t), d) == d ||
            throw(ArgumentError("cluster $c is treated at time $t for some observations and not for others"))
    end
    started = Dict{String, Int}()
    for (c, t, d) in zip(clusters, times, treated)
        d == 1 && (started[c] = min(get(started, c, t), t))
    end
    exposures = [d == 1 ? t - started[c] + 1 : 0 for (c, t, d) in zip(clusters, times, treated)]
    levels = sort!(unique(filter(>(0), exposures)))
    isempty(levels) && throw(ArgumentError("no observation is treated"))
    periods = sort!(unique(times))

    y = Float64[column(outcome)[k] for k in rows]
    weights = trials === nothing ? nothing : Float64[column(trials)[k] for k in rows]
    if weights !== nothing
        any(y .> weights) && throw(ArgumentError("the successes are at most the trials"))
        y ./= weights
    elseif family == "binomial"
        all(v -> v in (0, 1), y) || throw(ArgumentError("a binomial outcome without trials is 0 or 1"))
    end

    data = Dict{Symbol, AbstractVector}(:y => y, :cluster => clusters, :ij => string.(clusters, "-", times),
        :exposure => string.(exposures))
    individual === nothing || (data[:individual] = [string(column(individual)[k]) for k in rows])
    fixed = Symbol[]
    for t in periods
        name = Symbol("time_", t)
        data[name] = Float64.(times .== t)
        push!(fixed, name)
    end
    if effect == "ETI"
        for e in levels
            name = Symbol("exposure_", e)
            data[name] = Float64.(exposures .== e)
            push!(fixed, name)
        end
    else
        data[:treatment] = Float64.(treated)
        push!(fixed, :treatment)
    end
    groups = individual === nothing ? [:cluster, :ij] : [:cluster, :ij, :individual]
    fit = mixed_fit((; data...), fixed, groups, effect == "TEH" ? :treatment : nothing; family, weights)

    span = estimand == "PTE" ? (exposure:exposure) : exposure === nothing ? (levels[1]:levels[end]) : exposure
    chosen = filter(in(span), levels)
    isempty(chosen) && throw(ArgumentError("no observation is at exposure time $(estimand == "PTE" ? exposure : span)"))
    position(name) = findfirst(==(name), fit.names)
    if effect == "ETI"
        at = [position(Symbol("exposure_", e)) for e in levels]
        by_exposure = (; exposure = levels, effect = fit.coef[at], se = sqrt.(max.(diag(fit.vcov)[at], 0)))
        m = [e in chosen ? 1 / length(chosen) : 0.0 for e in levels]
        estimate = dot(m, fit.coef[at])
        se = sqrt(max(m' * fit.vcov[at, at] * m, 0))
    else
        j = position(:treatment)
        average, variance = fit.coef[j], fit.vcov[j, j]
        if effect == "TEH"
            modes = [fit.modes[string(e)] for e in levels]
            by_exposure = (; exposure = levels, effect = average .+ first.(modes), se = sqrt.(variance .+ last.(modes)))
        else
            by_exposure = (; exposure = levels, effect = fill(average, length(levels)), se = fill(sqrt(variance), length(levels)))
        end
        if effect == "TEH" && estimand == "PTE"
            k = findfirst(==(exposure), levels)
            estimate, se = by_exposure.effect[k], by_exposure.se[k]
        else
            estimate, se = average, sqrt(variance)
        end
    end
    sd = fit.sd
    icc = family == "gaussian" ?
        get(sd, :cluster, 0.0)^2 / (get(sd, :cluster, 0.0)^2 + get(sd, :ij, 0.0)^2 + get(sd, :individual, 0.0)^2 + fit.residual^2) :
        nothing
    SteppedWedgeResult(effect, estimand, estimand == "PTE" ? exposure : chosen[1]:chosen[end], estimate, se,
        (estimate - 1.96 * se, estimate + 1.96 * se), by_exposure,
        (; cluster = sd[:cluster], cluster_time = sd[:ij], individual = get(sd, :individual, nothing),
            exposure = get(sd, :exposure, nothing), residual = fit.residual),
        icc, (; observations = length(rows), clusters = length(unique(clusters)), times = length(periods),
            sequences = length(unique(values(started))), family), fit.model)
end

const EFFECT_NAMES = Dict("IT" => "immediate treatment", "ETI" => "exposure time indicators",
    "TEH" => "treatment effect heterogeneity")

function report(io::IO, result::SteppedWedgeResult)
    s, sd = result.sample, result.sd
    println(io, "\nStepped wedge: ", s.observations, " observations of ", s.clusters, " clusters over ", s.times, " times, ",
        s.sequences, " sequences; ", s.family, "; ", EFFECT_NAMES[result.effect], " (", result.effect, ")")
    println(io, "Standard deviations: clusters ", format(sd.cluster), ", clusters at a time ", format(sd.cluster_time),
        sd.individual === nothing ? "" : ", people $(format(sd.individual))",
        sd.exposure === nothing ? "" : ", effects by exposure time $(format(sd.exposure))",
        sd.residual === nothing ? "" : ", residual $(format(sd.residual))",
        result.icc === nothing ? "" : "; intracluster correlation $(format(result.icc))")
    b = result.by_exposure
    @printf(io, "%10s   %9s %8s\n", "exposure", "effect", "se")
    for j in eachindex(b.exposure)
        @printf(io, "%10d   %9s %8s\n", b.exposure[j], format(b.effect[j]), format(b.se[j]))
    end
    times = result.exposure isa Integer ? "exposure time $(result.exposure)" :
        "exposure times $(first(result.exposure)) to $(last(result.exposure))"
    println(io, result.estimand, " (", times, "): ", format(result.estimate), ", se ", format(result.se), ", 95% interval ",
        format(result.interval[1]), " to ", format(result.interval[2]))
end

"""
    write_results(result::SteppedWedgeResult, folder)

Writes the effects by exposure time (`sw-effects.csv`) and the estimate with its standard error, 95% interval, the
standard deviations of the random effects and the intracluster correlation (`sw-estimate.csv`) to `folder`; returns
the files.
"""
function write_results(result::SteppedWedgeResult, folder::AbstractString)
    mkpath(folder)
    [CSV.write(joinpath(folder, name), table) for (name, table) in result_tables(result)]
end

function result_tables(result::SteppedWedgeResult)
    value(x) = something(x, missing)
    sd = result.sd
    estimate = (; effect = [result.effect], estimand = [result.estimand], exposure = [string(result.exposure)],
        estimate = [result.estimate], se = [result.se], low = [result.interval[1]], high = [result.interval[2]],
        sd_cluster = [sd.cluster], sd_cluster_time = [sd.cluster_time], sd_individual = [value(sd.individual)],
        sd_exposure = [value(sd.exposure)], sd_residual = [value(sd.residual)], icc = [value(result.icc)])
    ["sw-effects.csv" => result.by_exposure, "sw-estimate.csv" => estimate]
end

# A panel study's result by one of its estimators, and its results, one for each estimator
const PanelStudyResult = Union{PanelResult, CallawaySantAnnaResult, SteppedWedgeResult}

"""
    result_tables(result)

The tables of a panel study's result by one of its estimators, each with the name of the file `write_results` writes
it to: `panel-effects.csv`, `cs-overall.csv`, `sw-estimate.csv` and so on.
"""
result_tables

report(io::IO, results::Vector{PanelStudyResult}) = foreach(r -> report(io, r), results)
write_results(results::Vector{PanelStudyResult}, folder::AbstractString) =
    reduce(vcat, [write_results(r, folder) for r in results])
