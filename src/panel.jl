# Panel event studies: units first treated at different times, against units not yet or never treated

"""
A panel event study's folder: `study.toml` (`source = "panel"`) names the table (`data`, a CSV in the folder by
default `panel.csv`), its columns and the event times, as the keywords of `panel_event_study`, and the estimators
(`estimators`, "panel_event_study" and "callaway_santanna"; the first by default). A `[callaway_santanna]` table holds
the keywords of `callaway_santanna` that differ from the others, such as its control units.
"""
struct PanelStudy
    folder::String
    title::String
    data::String
    estimators::Vector{String}
    options::Dict{Symbol, Any}
    callaway_options::Dict{Symbol, Any}
end

const PANEL_OPTIONS = (:unit, :time, :outcome, :cohort, :treatment, :cluster, :weights, :covariates, :leads, :lags,
    :reference, :control, :small_sample, :fixed_shares, :distribution, :uniform, :multiplier, :wild, :wild_weights,
    :draws, :seed)

const CALLAWAY_OPTIONS = (:unit, :time, :outcome, :cohort, :treatment, :cluster, :weights, :covariates, :leads, :lags,
    :control, :distribution, :base_period, :method, :anticipation, :balance, :unbalanced, :bootstrap, :draws, :uniform, :seed)

const PANEL_ESTIMATORS = ("panel_event_study", "callaway_santanna")

function read_panel_study(folder, settings)
    data = get(settings, "data", "panel.csv")
    estimators = String.(get(settings, "estimators", ["panel_event_study"]))
    for e in estimators
        e in PANEL_ESTIMATORS || throw(ArgumentError("estimators: \"$e\" is not one of $(join(repr.(PANEL_ESTIMATORS), ", "))"))
    end
    isempty(estimators) && throw(ArgumentError("estimators: none given"))
    options = Dict{Symbol, Any}(Symbol(key) => value for (key, value) in settings if Symbol(key) in PANEL_OPTIONS)
    callaway = Dict{Symbol, Any}(Symbol(key) => value for (key, value) in settings if Symbol(key) in CALLAWAY_OPTIONS)
    for (key, value) in get(settings, "callaway_santanna", Dict{String, Any}())
        Symbol(key) in CALLAWAY_OPTIONS || throw(ArgumentError("[callaway_santanna]: \"$key\" is not a keyword of callaway_santanna"))
        callaway[Symbol(key)] = value
    end
    PanelStudy(folder, get(settings, "title", basename(folder)), isabspath(data) ? data : joinpath(folder, data),
        unique(estimators), options, callaway)
end

"""
    analyse(study::PanelStudy)

The panel event study of the table in the study's folder by each of its estimators, as `panel_event_study` and
`callaway_santanna` with the keywords of its `study.toml`: a result each, in the order of `estimators`.
"""
function analyse(study::PanelStudy; now = nothing)
    table = CSV.File(study.data)
    Union{PanelResult, CallawaySantAnnaResult}[e == "panel_event_study" ? panel_event_study(table; study.options...) :
        callaway_santanna(table; study.callaway_options...) for e in study.estimators]
end

"""A panel event study's results: the sample, the effects by event time, each cohort's, and the tests of pre-trends."""
struct PanelResult
    # Observations, units, clusters, the cohorts, units never treated, the control units and the units dropped as
    # treated from the first time
    sample::NamedTuple
    # An event time a row: the observations of treated units at it, and for two-way fixed effects (twfe) and Sun and
    # Abraham (sa) the effect, its standard error, t, p and 95% band
    effects::NamedTuple
    # A cohort and an event time a row: Sun and Abraham's effect of the cohort at the event time, its standard error and
    # its weight in the event time's effect
    cohorts::NamedTuple
    # An estimator a row: the F test that the effects before treatment are all zero
    pretrends::NamedTuple
end

# A time or a cohort: a whole number of periods (a year, a quarter or a month numbered)
function whole(x)
    x isa Integer && return Int(x)
    x isa Real && isinteger(x) && return Int(x)
    throw(ArgumentError("times and cohorts are whole numbers of periods (years, or quarters or months numbered), not $x"))
end

"""
The time each unit is first treated, by unit, or nothing for a unit never treated: from a column of that time
(missing or 0 for never), or from the first time a treatment column is not zero.
"""
function first_treated(units, times, cohorts, treatments)
    first = Dict{String, Union{Int, Nothing}}()
    if cohorts !== nothing
        for (unit, cohort) in zip(units, cohorts)
            ismissing(unit) && continue
            e = ismissing(cohort) || cohort == 0 ? nothing : whole(cohort)
            key = string(unit)
            haskey(first, key) && first[key] != e && throw(ArgumentError("unit $key has two cohorts, $(first[key]) and $e"))
            first[key] = e
        end
    else
        for (unit, time, treated) in zip(units, times, treatments)
            (ismissing(unit) || ismissing(time)) && continue
            key = string(unit)
            e = get(first, key, nothing)
            first[key] = !ismissing(treated) && treated != 0 ? (e === nothing ? whole(time) : min(e, whole(time))) : e
        end
    end
    first
end

"""
Columns `X` and `y` less their weighted means by unit and by time, in turn until neither moves (the two-way fixed
effects, absorbed by alternating projections).
"""
function absorb!(M, units, times, w)
    groups = ((units, maximum(units)), (times, maximum(times)))
    totals = [accumulate_by(index, count, w) for (index, count) in groups]
    for column in eachcol(M)
        scale = max(1.0, maximum(abs, column; init = 0.0))
        for _ in 1:100_000
            change = 0.0
            for ((index, count), total) in zip(groups, totals)
                sums = accumulate_by(index, count, w .* column)
                for k in eachindex(column)
                    m = sums[index[k]] / total[index[k]]
                    column[k] -= m
                    change = max(change, abs(m))
                end
            end
            change < 1e-13 * scale && break
        end
    end
    M
end

function accumulate_by(index, count, values)
    sums = zeros(count)
    for k in eachindex(index)
        sums[index[k]] += values[k]
    end
    sums
end

"""
Weighted least squares of `y` on `X` with unit and time fixed effects, and the variance clustered by `clusters`, with
Stata's small-sample factor G / (G - 1) (N - 1) / (N - K) unless `small_sample` is false; K counts the time effects, the constant and, unless each
unit is within a cluster, the unit effects. Columns that are constant within the fixed effects or a combination of
others are left out: their coefficients are NaN. Also each cluster's influence on the coefficients, a row a cluster,
whose cross-product times the small-sample factor is V, for the multiplier bootstrap.
"""
function fixed_effects_regression(y, X, units, times, clusters, w, small_sample)
    M = absorb!(hcat(y, X), units, times, w)
    yd, Xd = M[:, 1], M[:, 2:end]
    p = size(X, 2)
    b, V = fill(NaN, p), fill(NaN, p, p)
    G = maximum(clusters)
    if p > 0
        factor = qr(sqrt.(w) .* Xd, ColumnNorm())
        R = abs.(diag(factor.R))
        rank = count(R .> maximum(R; init = 0.0) * 1e-9)
        kept = sort(factor.p[1:rank])
    else
        kept = Int[]
    end
    Xk = Xd[:, kept]
    A = Symmetric(Xk' * (w .* Xk))
    bk = A \ (Xk' * (w .* yd))
    u = yd - Xk * bk
    bread = inv(A)
    scores = zeros(G, length(kept))
    for k in eachindex(u)
        scores[clusters[k], :] .+= w[k] * u[k] .* Xk[k, :]
    end
    nested = all(length(unique(clusters[units .== i])) == 1 for i in unique(units))
    N = length(y)
    K = length(kept) + maximum(times) - 1 + 1 + (nested ? 0 : maximum(units) - 1)
    b[kept] = bk
    influence = fill(NaN, G, p)
    influence[:, kept] = scores * bread
    V[kept, kept] = (small_sample ? G / (G - 1) * (N - 1) / (N - K) : 1.0) * (influence[:, kept]' * influence[:, kept])
    (; b, V, clusters = G, influence)
end

"""
    panel_event_study(table; unit, time, outcome, cohort = nothing, treatment = nothing, cluster = unit,
        weights = nothing, covariates = String[], leads = 5, lags = 5, reference = -1, control = "auto", small_sample = true,
        fixed_shares = false, distribution = "t", uniform = false, multiplier = false, wild = false,
        wild_weights = "rademacher", draws = 10000, seed = nothing)

The effects of a treatment that units of a panel start at different times, by event time (the time since the unit was
first treated), from `-leads` to `lags`, against the `reference` event time:

- two-way fixed effects (TWFE): the outcome on unit and time effects and an indicator of each event time;
- Sun and Abraham (2021): the effect of each cohort (the units first treated at the same time) at each event time,
  against the control units, averaged over the cohorts with the shares of their observations at the event time.

With effects that differ from cohort to cohort, those of TWFE mix the effects of other event times in, with weights
that can be negative; Sun and Abraham's do not.

`table` is any Tables.jl table, a row a unit and a time; `unit`, `time` and `outcome` name its columns. Times are whole
numbers: years, or quarters or months numbered. When each unit was first treated comes from `cohort`, a column of that
time (missing or 0 for a unit never treated), or from `treatment`, the first time its column is not zero. Event times
before `-leads` count as `-leads` and those after `lags` as `lags`. Units treated from the first time in the table are
left out: they have no time before treatment.

`control` is "never" (the units never treated), "last" (the units treated last, with the times from their treatment
left out), or "auto": the units never treated when there are any, else those treated last.

Standard errors are clustered by `cluster` (by unit by default), with Stata's small-sample factor, as R's fixest,
unless `small_sample` is false, as R's did. With `distribution = "t"`, p-values and the 95% band come from Student's t
with the clusters less one degrees of freedom, as Stata and fixest; with `distribution = "normal"`, from the normal, as
R's did, which is right as the clusters grow in number. Sun and Abraham's variance adds that of the cohorts' shares, clustered as
the cohorts' effects; Stata's eventstudyinteract takes the shares' variance observation by observation, which gives
the same variance at each event time when each unit is a cluster and no event time is beyond the ends, but no
covariance between event times; R's did takes it as here. With `fixed_shares`, the shares are held fixed, as R's
fixest. To compare with R: fixest's effects and standard errors with `fixed_shares = true`, did's with
`small_sample = false` and `distribution = "normal"`. The pre-trends test is the Wald test that the effects before the
reference event time are all zero: F with the clusters less one degrees of freedom under t, chi-squared under the normal. `weights` weighs the
observations; `covariates` are columns added to both regressions.

With `uniform`, the 95% bands of each estimator's effects are simultaneous: they cover all the event times together
95% of the time, as an event-study plot needs. The critical value is the 95% quantile of the largest |t| over the event
times, in `draws` draws from the effects' estimated joint distribution (sup-t bands, Montiel Olea and Plagborg-Moller
2019): normal, or under `distribution = "t"` Student's t with the clusters less one degrees of freedom, so that with one
event time it is the pointwise value. With `multiplier`, the critical value comes from R's did's multiplier bootstrap
instead: each cluster's influence on the effects weighed by +1 or -1 in each draw, the largest |t| over the event times
each over its draws' interquartile range, and the pointwise value when it is smaller, as did's `aggte`. `seed` makes the
draws repeatable.

With `wild`, after `using WildBootTests`, each effect also has the p-value of the wild cluster bootstrap (Roodman,
MacKinnon, Nielsen and Webb 2019), with the null imposed, and its 95% interval by inverting the test, and the
pre-trends test its own p-value, for few clusters, where the clustered standard errors are unreliable. Each draw weighs
each cluster's residuals by +1 or -1 (`wild_weights = "rademacher"`; every pattern of signs when there are no more than
`draws`) or by Webb's six values (`"webb"`, better with fewer than about 12 clusters). Sun and Abraham's hold the
cohorts' shares fixed. The bootstrap is that of R's fwildclusterboot and Stata's boottest.
"""
function panel_event_study(table; unit, time, outcome, cohort = nothing, treatment = nothing, cluster = unit,
        weights = nothing, covariates = String[], leads::Integer = 5, lags::Integer = 5, reference::Integer = -1,
        control::AbstractString = "auto", small_sample::Bool = true, fixed_shares::Bool = false,
        distribution::AbstractString = "t", uniform::Bool = false, multiplier::Bool = false, wild::Bool = false,
        wild_weights::AbstractString = "rademacher", draws::Integer = 10000, seed = nothing)
    draws >= 100 || throw(ArgumentError("draws is at least 100"))
    multiplier && !uniform && throw(ArgumentError("multiplier gives the critical value of simultaneous bands: take uniform = true"))
    wild_weights in ("rademacher", "webb") || throw(ArgumentError("wild_weights is \"rademacher\" or \"webb\", not \"$wild_weights\""))
    wild && Base.get_extension(@__MODULE__, :PollisEventStudiesWildBootTestsExt) === nothing &&
        throw(ArgumentError("the wild cluster bootstrap needs WildBootTests.jl: using WildBootTests"))
    distribution in ("t", "normal") || throw(ArgumentError("distribution is \"t\" or \"normal\", not \"$distribution\""))
    (cohort === nothing) == (treatment === nothing) &&
        throw(ArgumentError("name the column of the time each unit is first treated (cohort) or of its treatment (treatment), one of the two"))
    -leads <= reference <= lags || throw(ArgumentError("the reference event time, $reference, is not within -$leads to $lags"))
    control in ("auto", "never", "last") || throw(ArgumentError("control is \"never\", \"last\" or \"auto\", not \"$control\""))
    columns = Tables.columns(table)
    column(name) = Tables.getcolumn(columns, Symbol(name))
    units, times = column(unit), column(time)
    treated_from = first_treated(units, times, cohort === nothing ? nothing : column(cohort),
        treatment === nothing ? nothing : column(treatment))

    # The observations with all that is needed
    needed = [column(name) for name in (unit, time, outcome, cluster, covariates...)]
    weighed = weights === nothing ? nothing : column(weights)
    rows = [k for k in eachindex(units) if all(!ismissing(c[k]) for c in needed) &&
        (weighed === nothing || (!ismissing(weighed[k]) && weighed[k] > 0))]
    isempty(rows) && throw(ArgumentError("no observation has the unit, time, outcome, cluster and covariates"))
    start = minimum(whole(times[k]) for k in rows)
    cohort_of(k) = treated_from[string(units[k])]
    always = Set(string(units[k]) for k in rows if cohort_of(k) !== nothing && cohort_of(k) <= start)
    filter!(k -> !(string(units[k]) in always), rows)
    never = any(cohort_of(k) === nothing for k in rows)
    control == "auto" && (control = never ? "never" : "last")
    control == "never" && !never && throw(ArgumentError("no unit is never treated: take control = \"last\""))
    cohorts = sort!(unique(cohort_of(k) for k in rows if cohort_of(k) !== nothing))
    control_cohort = control == "last" ? (isempty(cohorts) ? nothing : cohorts[end]) : nothing
    # Against the units treated last, the times from their treatment have no control
    control_cohort === nothing || filter!(k -> whole(times[k]) < control_cohort, rows)
    treated = filter(e -> e != control_cohort, cohorts)
    isempty(treated) && throw(ArgumentError("no unit is treated within the table, apart from the control units"))

    index(values) = (levels = sort!(unique(values)); positions = Dict(v => i for (i, v) in enumerate(levels)); [positions[v] for v in values])
    unit_index = index([string(units[k]) for k in rows])
    time_index = index([whole(times[k]) for k in rows])
    cluster_index = index([string(column(cluster)[k]) for k in rows])
    w = weighed === nothing ? ones(length(rows)) : Float64[weighed[k] for k in rows]
    y = Float64[column(outcome)[k] for k in rows]
    Z = Float64[c[k] for k in rows, c in [column(name) for name in covariates]]
    e = [cohort_of(k) for k in rows]
    event = [x === nothing ? nothing : clamp(whole(times[k]) - x, -leads, lags) for (k, x) in zip(rows, e)]
    event_times = [l for l in -leads:lags if l != reference]

    # Two-way fixed effects: an indicator of each event time, for every treated unit
    X = Float64[event[k] == l for k in eachindex(rows), l in event_times]
    twfe = fixed_effects_regression(y, hcat(X, Z), unit_index, time_index, cluster_index, w, small_sample)

    # Sun and Abraham: an indicator of each cohort at each event time, but for the control units
    pairs = [(c, l) for c in treated for l in event_times if any(e[k] == c && event[k] == l for k in eachindex(rows))]
    D = Float64[e[k] == c && event[k] == l for k in eachindex(rows), (c, l) in pairs]
    sa = fixed_effects_regression(y, hcat(D, Z), unit_index, time_index, cluster_index, w, small_sample)
    delta = sa.b[1:length(pairs)]
    V_delta = sa.V[1:length(pairs), 1:length(pairs)]

    # Each event time's effect: the cohorts' at it, weighed by their shares of its observations, and its variance:
    # that of the cohorts' effects and that of their shares, both clustered
    L = length(event_times)
    shares = zeros(length(pairs), L)
    influence = zeros(sa.clusters, L)
    observations = zeros(Int, L)
    for (j, l) in enumerate(event_times)
        at = [k for k in eachindex(rows) if e[k] in treated && event[k] == l]
        observations[j] = length(at)
        isempty(at) && continue
        total = sum(w[at])
        for (i, (c, m)) in enumerate(pairs)
            m == l && (shares[i, j] = sum(w[k] for k in at if e[k] == c; init = 0.0) / total)
        end
        # The influence of each cluster on the shares, through the cohorts' effects at the event time
        for k in at
            for (i, (c, m)) in enumerate(pairs)
                m == l || continue
                influence[cluster_index[k], j] += w[k] * ((e[k] == c) - shares[i, j]) * delta[i] / total
            end
        end
    end
    used(j) = [i for i in eachindex(pairs) if shares[i, j] > 0]
    beta = [isempty(used(j)) || any(isnan, delta[used(j)]) ? NaN : sum(shares[i, j] * delta[i] for i in used(j)) for j in 1:L]
    V_beta = shares' * replace(V_delta, NaN => 0.0) * shares + (fixed_shares ? 0.0 : 1.0) * (influence' * influence)
    for j in 1:L
        isnan(beta[j]) && (V_beta[j, :] .= NaN; V_beta[:, j] .= NaN)
    end

    df = sa.clusters - 1
    normal = distribution == "normal"
    critical = quantile(normal ? Normal() : TDist(df), 0.975)
    probability(t) = ismissing(t) || !isfinite(t) ? missing : normal ? 2 * ccdf(Normal(), abs(t)) : pvalue(t, df)
    function estimates(b, V, c = critical)
        se = [isnan(V[j, j]) ? NaN : sqrt(V[j, j]) for j in eachindex(b)]
        t = b ./ se
        (; effect = nan_missing.(b), se = nan_missing.(se), t = nan_missing.(t), p = [probability(nan_missing(x)) for x in t],
            low = nan_missing.(b .- c .* se), high = nan_missing.(b .+ c .* se))
    end
    # Simultaneous bands: the sup-t from the effects' covariance, or did's multiplier bootstrap on each cluster's
    # influence on them (Sun and Abraham's: the cohorts' effects' through the shares, and the shares' own)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    function simultaneous(V, Psi)
        uniform || return critical
        if multiplier
            known = [j for j in 1:L if !isnan(V[j, j])]
            c = multiplier_critical(multiplier_draws(rng, Psi[:, known], draws))
            return c === nothing || c < quantile(Normal(), 0.975) ? critical : c
        end
        c = supt_critical(rng, V, draws, normal ? nothing : df)
        c === nothing ? critical : c
    end
    sa_influence = replace(sa.influence[:, 1:length(pairs)], NaN => 0.0) * shares + (fixed_shares ? 0.0 : 1.0) * influence
    twfe_critical = simultaneous(twfe.V[1:L, 1:L], twfe.influence[:, 1:L])
    sa_critical = simultaneous(V_beta, sa_influence)
    first_estimates = estimates(twfe.b[1:L], twfe.V[1:L, 1:L], twfe_critical)
    second_estimates = estimates(beta, V_beta, sa_critical)
    effects = (; event_time = event_times, observations,
        (Symbol("twfe_", k) => v for (k, v) in pairs_of(first_estimates))...,
        (Symbol("sa_", k) => v for (k, v) in pairs_of(second_estimates))...)
    cohort_se = [isnan(V_delta[i, i]) ? NaN : sqrt(V_delta[i, i]) for i in eachindex(pairs)]
    cohort_table = (; cohort = first.(pairs), event_time = last.(pairs), effect = nan_missing.(delta), se = nan_missing.(cohort_se),
        weight = [shares[i, findfirst(==(l), event_times)] for (i, (_, l)) in enumerate(pairs)])
    before = [j for j in 1:L if event_times[j] < reference]
    # The Wald test: F(q, clusters less one) under t, chi-squared(q), q times F, under the normal
    tests = [wald(twfe.b[1:L], twfe.V[1:L, 1:L], before), wald(beta, V_beta, before)]
    statistic = [ismissing(F) ? missing : normal ? q * F : F for (F, q) in tests]
    pretrends = (; estimator = ["Two-way fixed effects", "Sun and Abraham"], test = fill(normal ? "chi2" : "F", 2),
        statistic, df1 = last.(tests), df2 = fill(normal ? missing : df, 2),
        p = [ismissing(x) ? missing : normal ? ccdf(Chisq(q), x) : ccdf(FDist(q, df), x) for (x, (_, q)) in zip(statistic, tests)])
    # The wild cluster bootstrap: the regression's columns with a coefficient, the time effects as indicators, the unit
    # effects absorbed; each effect a combination of the columns (Sun and Abraham's, the shares held fixed)
    wild_draws = 0
    if wild
        time_dummies = Float64[time_index[k] == t for k in eachindex(rows), t in 2:maximum(time_index)]
        function wild_tests(regression, regressors, combination)
            kept = findall(!isnan, regression.b)
            design = hcat(regressors[:, kept], time_dummies)
            testable(j) = all(iszero, combination(j)[setdiff(axes(regressors, 2), kept)]) && any(!iszero, combination(j))
            restriction(js) = hcat(reduce(vcat, [combination(j)[kept]' for j in js]), zeros(length(js), size(time_dummies, 2)))
            test(js, interval) = wild_test(restriction(js), y, design, unit_index, cluster_index, w; draws, rng,
                auxiliary = wild_weights, interval)
            single = [testable(j) ? test([j], true) : nothing for j in 1:L]
            joint = filter(testable, before)
            pre = isempty(joint) ? missing : test(joint, false).p
            wild_draws = maximum((t.draws for t in single if t !== nothing); init = wild_draws)
            ([t === nothing ? missing : t.p for t in single], [t === nothing ? missing : t.interval[1] for t in single],
                [t === nothing ? missing : t.interval[2] for t in single], pre)
        end
        unit_vector(j) = (c = zeros(size(X, 2) + size(Z, 2)); c[j] = 1.0; c)
        twfe_wild = wild_tests(twfe, hcat(X, Z), unit_vector)
        sa_wild = wild_tests(sa, hcat(D, Z), j -> vcat(isnan(beta[j]) ? zeros(length(pairs)) : shares[:, j], zeros(size(Z, 2))))
        effects = merge(effects, (; twfe_wild_p = twfe_wild[1], twfe_wild_low = twfe_wild[2], twfe_wild_high = twfe_wild[3],
            sa_wild_p = sa_wild[1], sa_wild_low = sa_wild[2], sa_wild_high = sa_wild[3]))
        pretrends = merge(pretrends, (; wild_p = [twfe_wild[4], sa_wild[4]]))
    end
    sample = (; observations = length(rows), units = maximum(unit_index), clusters = sa.clusters, cohorts = treated,
        never_treated = length(unique(unit_index[k] for k in eachindex(rows) if e[k] === nothing)),
        control = control == "never" ? "never treated" : "treated last ($control_cohort), before their treatment",
        dropped = length(always), leads, lags, reference, uniform, multiplier, draws = uniform ? draws : 0, wild,
        wild_weights, wild_draws,
        critical = (; twfe = twfe_critical, sa = sa_critical))
    PanelResult(sample, effects, cohort_table, pretrends)
end

nan_missing(x) = isnan(x) ? missing : x
pairs_of(named) = zip(keys(named), values(named))

"""The F statistic that the effects `which` are all zero, and its degrees of freedom; missing without any."""
function wald(b, V, which)
    which = [j for j in which if !isnan(b[j])]
    isempty(which) && return (missing, 0)
    v = b[which]
    (v' * (Symmetric(V[which, which]) \ v) / length(which), length(which))
end

function report(io::IO, result::PanelResult)
    s = result.sample
    println(io, "\n", s.observations, " observations of ", s.units, " units in ", s.clusters, " clusters; cohorts ",
        join(s.cohorts, ", "), "; ", s.never_treated, " never treated; control: ", s.control,
        s.dropped > 0 ? "; $(s.dropped) units treated from the first time left out" : "")
    println(io, "Event times ", -s.leads, " to ", s.lags, " against ", s.reference, " (the ends take the times beyond them)")
    s.uniform && println(io, "95% bands simultaneous over the event times (", s.multiplier ? "did's multiplier bootstrap" : "sup-t",
        ", ", s.draws, " draws): critical values ",
        format(s.critical.twfe), " (TWFE), ", format(s.critical.sa), " (Sun and Abraham)")
    effects = result.effects
    @printf(io, "%6s %6s   %9s %8s %8s %7s   %9s %8s %8s %7s\n", "event", "obs", "TWFE", "se", "t", "p", "Sun-Abr.", "se", "t", "p")
    for j in eachindex(effects.event_time)
        @printf(io, "%6d %6d   %9s %8s %8s %7s   %9s %8s %8s %7s\n", effects.event_time[j], effects.observations[j],
            format(effects.twfe_effect[j]), format(effects.twfe_se[j]), format(effects.twfe_t[j]), format(effects.twfe_p[j]),
            format(effects.sa_effect[j]), format(effects.sa_se[j]), format(effects.sa_t[j]), format(effects.sa_p[j]))
    end
    for k in eachindex(result.pretrends.estimator)
        p = result.pretrends
        degrees = p.test[k] == "F" ? "F($(p.df1[k]), $(p.df2[k]))" : "chi2($(p.df1[k]))"
        println(io, "Pre-trends, ", p.estimator[k], ": ", degrees, " = ", format(p.statistic[k]), ", p = ", format(p.p[k]),
            s.wild ? string("; wild bootstrap p = ", format(p.wild_p[k])) : "")
    end
    s.wild || return
    println(io, "Wild cluster bootstrap (", s.wild_weights == "webb" ? "Webb's" : "Rademacher", " weights, ", s.wild_draws,
        " draws by ", s.clusters, " clusters): p-values, the null imposed, and 95% intervals by inverting the test")
    @printf(io, "%6s   %7s %9s %9s   %7s %9s %9s\n", "event", "TWFE p", "low", "high", "S-A p", "low", "high")
    for j in eachindex(effects.event_time)
        @printf(io, "%6d   %7s %9s %9s   %7s %9s %9s\n", effects.event_time[j], format(effects.twfe_wild_p[j]),
            format(effects.twfe_wild_low[j]), format(effects.twfe_wild_high[j]), format(effects.sa_wild_p[j]),
            format(effects.sa_wild_low[j]), format(effects.sa_wild_high[j]))
    end
end

"""
    write_results(result::PanelResult, folder)

Writes the effects by event time (`panel-effects.csv`), each cohort's (`panel-cohorts.csv`) and the tests of
pre-trends (`panel-pretrends.csv`) to `folder`; returns the files.
"""
function write_results(result::PanelResult, folder::AbstractString)
    mkpath(folder)
    [CSV.write(joinpath(folder, name), table) for (name, table) in
        (("panel-effects.csv", result.effects), ("panel-cohorts.csv", result.cohorts), ("panel-pretrends.csv", result.pretrends))]
end

heading(study::PanelStudy) = string(study.title, ": panel event study of ", basename(study.data))
