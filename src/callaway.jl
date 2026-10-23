# Callaway and Sant'Anna (2021): the effect of each cohort at each time, from a difference in differences against
# the units never or not yet treated, aggregated by event time, by cohort and by calendar time. A port of R's did
# (att_gt and aggte) and of DRDID's 2x2 estimators for panels (drdid_panel, std_ipw_did_panel, reg_did_panel) and, for
# unbalanced panels, for repeated cross sections (drdid_rc, std_ipw_did_rc, reg_did_rc).

"""
Callaway and Sant'Anna's results: the sample, each cohort's effect at each time (`cells`), their aggregations by event
time (`effects`), by cohort (`cohorts`) and by calendar time (`calendar`), the overall averages (`overall`) and the
test of pre-trends.
"""
struct CallawaySantAnnaResult
    # Units, clusters, times, cohorts, units never treated, the control units, the base period, the method, the units
    # dropped (treated from the first time), the units not observed at every time (left out, or kept when panel is
    # "unbalanced")
    sample::NamedTuple
    # A cohort and a time a row: ATT(g, t), its standard error, t, p and 95% band; post is whether the cohort is
    # treated at the time
    cells::NamedTuple
    # An event time a row: the cohorts' effects at it, weighed by the cohorts' sizes
    effects::NamedTuple
    # A cohort a row: the mean of its effects from its treatment on
    cohorts::NamedTuple
    # A time a row: the effects of the cohorts treated by then, weighed by their sizes
    calendar::NamedTuple
    # An aggregation a row: the effects from treatment on, averaged over the cells (simple), the event times, the
    # cohorts and the times
    overall::NamedTuple
    # The Wald test that the effects of the cells before treatment are all zero
    pretrends::NamedTuple
end

"""
    callaway_santanna(table; unit, time, outcome, cohort = nothing, treatment = nothing, cluster = unit,
        weights = nothing, covariates = String[], control = "auto", base_period = "varying", method = "dr",
        anticipation = 0, leads = nothing, lags = nothing, balance = nothing, unbalanced = false, distribution = "normal",
        bootstrap = false, draws = 1000, uniform = bootstrap, seed = nothing)

Callaway and Sant'Anna's (2021) effects of a treatment that units of a panel start at different times. For each cohort
g (the units first treated at time g) and each time t, ATT(g, t) is a difference in differences: the change in the
outcome of the cohort from the base period to t, less that of the control units.

`table`, `unit`, `time`, `outcome`, `cohort`, `treatment`, `cluster` and `weights` are as for `panel_event_study`; a
cohort after the table's last time counts as never treated. Units treated from the first time are left out, as are
units not observed at every time unless `unbalanced`.

- `control`: "never" (the units never treated), "notyet" (those never treated and those not yet treated by t and by
  the base period), "last" (those never treated and those treated last, with the times from their treatment left out)
  or "auto": "never" when there are units never treated, else "last".
- `base_period`: "varying", the time before t for the cells before treatment and the last time before the cohort's
  treatment for those after; or "universal", the last time before the cohort's treatment for all, whose cell is zero.
- `method`, with `covariates` (taken at the earlier of t and the base period): "dr" doubly robust (Sant'Anna and Zhao
  2020), "ipw" inverse probability weighting or "reg" outcome regression. Without covariates the three are the
  difference in the mean changes.
- `anticipation`: the times before its treatment a cohort may already respond.
- `unbalanced`: keeps the units not observed at every time, as did's `allow_unbalanced_panel`. Each ATT(g, t) then
  compares the observations at t and at the base period as repeated cross sections (DRDID's drdid_rc, std_ipw_did_rc
  and reg_did_rc), each with its own covariates and weight, a unit's influence being the sum of its observations'; the
  cohorts' sizes weigh each unit by the mean of its weights. When every unit is observed at every time, nothing
  changes.

The cells are aggregated by event time (from `-leads` to `lags`, or all of them; with `balance` = e, only the cohorts
seen e times after their treatment), by cohort and by calendar time, each weighed by the cohorts' sizes, and their
standard errors count the estimation of the weights. `lags` also bounds the cells averaged by cohort and overall.

Standard errors come from the influence functions, clustered by `cluster`, without a small-sample factor. p-values and
95% bands come from the normal (`distribution = "normal"`, the default, as did, the authors' R package) or from
Student's t with the clusters less one degrees of freedom (`distribution = "t"`, wider with few clusters), and the test
of pre-trends is chi-squared or F; there is none, as in did, when the covariance of the cells before treatment is
singular (more of them than clusters). The estimates, standard errors, p-values and bands are those of R's did with
`bstrap = FALSE, cband = FALSE`.

With `bootstrap = true`, the standard errors come from did's multiplier bootstrap instead: `draws` times, each
cluster's influence is weighed by +1 or -1 at random, and the standard error is the interquartile range of the draws
over that of the normal. With `uniform` (as did's `cband`, and by default with the bootstrap) the 95% bands of the
cells, and of the effects by event time, by cohort and by time, are simultaneous: each table's critical value is the
95% quantile of the largest t over its rows in the draws, so that the band covers all of them together 95% of the
time. Pointwise bands are kept when that value is not found or is below the normal's, and for the overall averages.
`seed` makes the draws repeatable.
"""
function callaway_santanna(table; unit, time, outcome, cohort = nothing, treatment = nothing, cluster = unit,
        weights = nothing, covariates = String[], control::AbstractString = "auto",
        base_period::AbstractString = "varying", method::AbstractString = "dr", anticipation::Integer = 0,
        leads = nothing, lags = nothing, balance = nothing, unbalanced::Bool = false,
        distribution::AbstractString = "normal", bootstrap::Bool = false,
        draws::Integer = 1000, uniform::Bool = bootstrap, seed = nothing)
    (cohort === nothing) == (treatment === nothing) &&
        throw(ArgumentError("name the column of the time each unit is first treated (cohort) or of its treatment (treatment), one of the two"))
    control in ("auto", "never", "notyet", "last") ||
        throw(ArgumentError("control is \"never\", \"notyet\", \"last\" or \"auto\", not \"$control\""))
    base_period in ("varying", "universal") || throw(ArgumentError("base_period is \"varying\" or \"universal\", not \"$base_period\""))
    method in ("dr", "ipw", "reg") || throw(ArgumentError("method is \"dr\", \"ipw\" or \"reg\", not \"$method\""))
    distribution in ("t", "normal") || throw(ArgumentError("distribution is \"t\" or \"normal\", not \"$distribution\""))
    anticipation >= 0 || throw(ArgumentError("anticipation is not negative"))
    uniform && !bootstrap && throw(ArgumentError("simultaneous bands (uniform) need the bootstrap: bootstrap = true"))
    draws >= 100 || throw(ArgumentError("draws is at least 100"))
    columns = Tables.columns(table)
    column(name) = Tables.getcolumn(columns, Symbol(name))
    units, times = column(unit), column(time)
    treated_from = first_treated(units, times, cohort === nothing ? nothing : column(cohort),
        treatment === nothing ? nothing : column(treatment))

    # The observations with all that is needed
    needed = [column(name) for name in (unit, time, outcome, cluster, covariates...)]
    weighed = weights === nothing ? nothing : column(weights)
    rows = [k for k in eachindex(units) if all(!ismissing(c[k]) for c in needed) && (weighed === nothing || !ismissing(weighed[k]))]
    weighed === nothing || all(weighed[k] >= 0 for k in rows) || throw(ArgumentError("weights are not negative"))
    isempty(rows) && throw(ArgumentError("no observation has the unit, time, outcome, cluster and covariates"))
    # A cohort after the last time is never treated, as R's did
    last_time = maximum(whole(times[k]) for k in rows)
    cohort_of(k) = (e = treated_from[string(units[k])]; e === nothing || e > last_time + anticipation ? 0 : e)
    never = any(cohort_of(k) == 0 for k in rows)
    control == "auto" && (control = never ? "never" : "last")
    control == "never" && !never && throw(ArgumentError("no unit is never treated: take control = \"last\" or \"notyet\""))

    # Against the units treated last (or, not yet treated, without units never treated), the times from their
    # treatment have no control; with control = "last" they are the units never treated
    latest = maximum(cohort_of(k) for k in rows)
    as_never = control == "last" ? latest : -1
    if control == "last" || (control == "notyet" && !never)
        filter!(k -> whole(times[k]) < latest - anticipation, rows)
    end
    g_of(k) = (e = cohort_of(k); e == as_never ? 0 : e)
    start = minimum(whole(times[k]) for k in rows)
    always = Set(string(units[k]) for k in rows if 0 < g_of(k) <= start + anticipation)
    filter!(k -> !(string(units[k]) in always), rows)

    # A balanced panel, the units seen at every time, once; or, unbalanced when some are not, all the units
    time_list = sort!(unique(whole(times[k]) for k in rows))
    T = length(time_list)
    seen = Dict{String, Vector{Int}}()
    for k in rows
        push!(get!(seen, string(units[k]), Int[]), k)
    end
    for (key, ks) in seen
        length(unique(whole(times[k]) for k in ks)) == length(ks) || throw(ArgumentError("unit $key is seen more than once at a time"))
    end
    gaps = count(ks -> length(ks) != T, values(seen))
    repeated = unbalanced && gaps > 0
    unit_list = sort!([key for (key, ks) in seen if repeated || length(ks) == T])
    isempty(unit_list) && throw(ArgumentError("no unit is seen at every time: take unbalanced = true"))
    n = length(unit_list)
    position = Dict(t => j for (j, t) in enumerate(time_list))
    Y = fill(NaN, n, T)
    observed = falses(n, T)
    Wt = ones(n, T)
    X = [ones(n, T)]
    for _ in covariates
        push!(X, zeros(n, T))
    end
    G = zeros(Int, n)
    never_treated = 0
    clusters = Vector{String}(undef, n)
    for (i, key) in enumerate(unit_list)
        for k in seen[key]
            j = position[whole(times[k])]
            Y[i, j] = column(outcome)[k]
            observed[i, j] = true
            weighed === nothing || (Wt[i, j] = weighed[k])
            for (c, covariate) in enumerate(covariates)
                X[c + 1][i, j] = column(covariate)[k]
            end
        end
        G[i] = g_of(first(seen[key]))
        never_treated += cohort_of(first(seen[key])) == 0
        clusters[i] = string(column(cluster)[first(seen[key])])
        all(string(column(cluster)[k]) == clusters[i] for k in seen[key]) ||
            throw(ArgumentError("unit $key changes cluster: clusters are fixed over time"))
    end
    groups = sort!(unique(g for g in G if g > start + anticipation))
    control == "notyet" && !never && filter!(<(latest), groups)
    isempty(groups) && throw(ArgumentError("no unit is treated within the table, apart from the control units"))

    # Each cohort at each time: ATT(g, t) and its influence function, n by the cells
    universal = base_period == "universal"
    shift = universal ? 0 : 1
    cell_group, cell_time, cell_att, cell_post = Int[], Int[], Float64[], Bool[]
    influence = Vector{Float64}[]
    for g in groups
        before = findlast(t -> t + anticipation < g, time_list)
        before === nothing && continue
        for s in 1:(T - shift)
            now = s + shift
            # The base period: the last time before the cohort's treatment, or, varying, the time before t before it
            base = universal || g <= time_list[now] ? before : s
            if universal && base == now
                push!(cell_group, g); push!(cell_time, time_list[now]); push!(cell_att, 0.0); push!(cell_post, false)
                push!(influence, zeros(n))
                continue
            end
            threshold = time_list[max(s, base) + shift] + anticipation
            treated = G .== g
            controls = control == "notyet" ? [G[i] == 0 || (G[i] > threshold && G[i] != g) for i in 1:n] : G .== 0
            kept = findall(treated .| controls)
            earlier = min(now, base)
            # Unbalanced, each observation at the two times on its own, with its covariates and weight
            obs = repeated ? [(i, j) for i in kept for j in (base, now) if observed[i, j]] : Tuple{Int, Int}[]
            D = Float64[treated[i] for (i, _) in obs]
            post = Float64[j == now for (_, j) in obs]
            att, inf = try
                if !repeated
                    two_by_two(Y[kept, now] .- Y[kept, base], Float64.(treated[kept]),
                        reduce(hcat, [x[kept, earlier] for x in X]), Wt[kept, earlier], method)
                elseif all(any((D .== d) .& (post .== p)) for d in (0, 1), p in (0, 1))
                    two_by_two_rc([Y[i, j] for (i, j) in obs], post, D, [x[i, j] for (i, j) in obs, x in X],
                        [Wt[i, j] for (i, j) in obs], method)
                else
                    # No effect, as did, without treated and control units at both times
                    (NaN, fill(NaN, length(obs)))
                end
            catch err
                err isa Union{LinearAlgebra.SingularException, LinearAlgebra.PosDefException, LinearAlgebra.LAPACKException,
                    ArgumentError, DomainError} || rethrow()
                (NaN, fill(NaN, repeated ? length(obs) : length(kept)))
            end
            full = zeros(n)
            if repeated
                # A unit's influence is the sum of its observations'
                for (k, (i, _)) in enumerate(obs)
                    full[i] += (n / length(obs)) * inf[k]
                end
            else
                full[kept] = (n / length(kept)) .* inf
            end
            isnan(att) && (full .= NaN)
            push!(cell_group, g); push!(cell_time, time_list[now]); push!(cell_att, att); push!(cell_post, g <= time_list[now])
            push!(influence, full)
        end
    end
    isempty(cell_att) && throw(ArgumentError("no cohort has a time before its treatment"))
    Psi = reduce(hcat, influence)

    # Standard errors from the influence functions, clustered
    cluster_index = let levels = sort!(unique(clusters)); positions = Dict(c => j for (j, c) in enumerate(levels)); [positions[c] for c in clusters] end
    C = maximum(cluster_index)
    by_cluster(v) = accumulate_by(cluster_index, C, v)
    standard_error(v) = any(isnan, v) ? NaN : (se = sqrt(sum(abs2, by_cluster(v))) / n; se <= sqrt(eps()) * 10 ? NaN : se)
    df = C - 1
    normal = distribution == "normal"
    critical = quantile(normal ? Normal() : TDist(df), 0.975)
    probability(t) = ismissing(t) || !isfinite(t) ? missing : normal ? 2 * ccdf(Normal(), abs(t)) : pvalue(t, df)
    function estimates(b, se, c = critical)
        t = b ./ se
        (; effect = nan_missing.(b), se = nan_missing.(se), t = nan_missing.(t), p = [probability(nan_missing(x)) for x in t],
            low = nan_missing.(b .- c .* se), high = nan_missing.(b .+ c .* se))
    end

    # did's multiplier bootstrap (mboot) on each cluster's sum of the influence functions (the columns of M)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    drawn(M) = multiplier_draws(rng, reduce(hcat, [by_cluster(M[:, k]) for k in axes(M, 2)]), draws)
    function bootstrap_error(v)
        any(isnan, v) && return NaN
        b = vec(drawn(reshape(v, :, 1)))
        degenerate(b) && return NaN
        se = normal_spread(b) * sqrt(C) / n
        se <= TINY ? NaN : se
    end
    error_of(v) = bootstrap ? bootstrap_error(v) : standard_error(v)
    # A table's simultaneous critical value, or the pointwise one; as did's att_gt for the cells (floor = false), and
    # aggte for the aggregations
    function simultaneous(M; floor = true)
        uniform || return critical
        c = multiplier_critical(drawn(M))
        c === nothing || (floor && c < quantile(Normal(), 0.975)) ? critical : c
    end

    analytic = [standard_error(Psi[:, j]) for j in axes(Psi, 2)]
    cell_se = bootstrap ? [isnan(analytic[j]) ? NaN : bootstrap_error(Psi[:, j]) for j in axes(Psi, 2)] : analytic
    cell_critical = simultaneous(Psi; floor = false)
    cells = (; cohort = cell_group, time = cell_time, event_time = cell_time .- cell_group, post = cell_post,
        estimates(cell_att, cell_se, cell_critical)...)

    # The test of pre-trends: the cells before treatment with a standard error, from the influence functions
    pre = [j for j in eachindex(cell_att) if cell_group[j] > cell_time[j] && !isnan(analytic[j])]
    S = isempty(pre) ? zeros(C, 0) : reduce(hcat, [by_cluster(Psi[:, j]) for j in pre])
    V = Symmetric(S' * S) / n^2
    # None without cells before treatment, or, as did, when their covariance is singular (more cells than clusters)
    pretrends = if isempty(pre) || inv(cond(V, 1)) <= eps()
        (; test = normal ? "chi2" : "F", statistic = missing, df1 = length(pre), df2 = normal ? missing : df, p = missing)
    else
        W = cell_att[pre]' * (V \ cell_att[pre])
        q = length(pre)
        normal ? (; test = "chi2", statistic = W, df1 = q, df2 = missing, p = ccdf(Chisq(q), W)) :
            (; test = "F", statistic = W / q, df1 = q, df2 = df, p = ccdf(FDist(q, df), W / q))
    end

    # Aggregations, over the cells with an effect, weighed by the cohorts' shares of the units (weighted) and with the
    # influence of the shares' estimation
    usable = findall(!isnan, cell_att)
    # Each unit's weight: at the first time, or, unbalanced, the mean of its weights, as did
    w1 = repeated ? [sum(Wt[i, observed[i, :]]) / count(observed[i, :]) for i in 1:n] : Wt[:, 1]
    share = Dict(g => sum(w1[i] for i in 1:n if G[i] == g; init = 0.0) / n for g in groups)
    max_e = lags === nothing ? Inf : lags
    min_e = leads === nothing ? -Inf : -leads
    function aggregate(which, weighted = true)
        isempty(which) && return (NaN, fill(NaN, n))
        p = [share[cell_group[j]] for j in which]
        total = sum(p)
        effect = sum(cell_att[j] * p[k] / total for (k, j) in enumerate(which))
        psi = Psi[:, which] * (p ./ total)
        if weighted
            centered = [w1[i] * (G[i] == cell_group[j]) - share[cell_group[j]] for i in 1:n, j in which]
            psi .+= (centered ./ total .- sum(centered; dims = 2) * (p' ./ total^2)) * cell_att[which]
        end
        (effect, psi)
    end
    function averaged(parts)
        isempty(parts) && return (NaN, fill(NaN, n))
        (sum(first, parts) / length(parts), sum(last, parts) / length(parts))
    end
    summary(parts, c = critical) = estimates(first.(parts), [error_of(last(p)) for p in parts], c)
    band(parts) = isempty(parts) ? critical : simultaneous(reduce(hcat, last.(parts)))

    included = trues(length(cell_att))
    if balance !== nothing
        included = [time_list[end] - g >= balance for g in cell_group]
    end
    event_times = sort!(unique(cell_time[j] - cell_group[j] for j in usable if included[j]))
    balance === nothing || filter!(e -> e <= balance && e >= balance - time_list[end] + time_list[1], event_times)
    filter!(e -> min_e <= e <= max_e, event_times)
    by_event = [aggregate([j for j in usable if included[j] && cell_time[j] - cell_group[j] == e]) for e in event_times]
    event_critical = band(by_event)
    effects = (; event_time = event_times, summary(by_event, event_critical)...)

    after(j) = cell_group[j] <= cell_time[j] <= cell_group[j] + max_e
    cohort_list = [g for g in groups if any(cell_group[j] == g && after(j) for j in usable)]
    by_cohort = [aggregate([j for j in usable if cell_group[j] == g && after(j)], false) for g in cohort_list]
    cohort_critical = band(by_cohort)
    cohort_table = (; cohort = cohort_list, summary(by_cohort, cohort_critical)...)

    calendar_times = [t for t in time_list if t >= minimum(cell_group[usable]; init = typemax(Int)) &&
        any(cell_time[j] == t && cell_group[j] <= t for j in usable)]
    by_time = [aggregate([j for j in usable if cell_time[j] == t && cell_group[j] <= t]) for t in calendar_times]
    calendar_critical = band(by_time)
    calendar = (; time = calendar_times, summary(by_time, calendar_critical)...)

    # Overall: the cells from treatment on weighed by the cohorts' sizes; the mean of the event times from 0; the
    # cohorts' effects weighed by their sizes; the mean of the times
    simple = aggregate([j for j in usable if after(j)])
    over_events = averaged([by_event[k] for k in eachindex(event_times) if event_times[k] >= 0])
    over_cohorts = if isempty(cohort_list)
        (NaN, fill(NaN, n))
    else
        p = [share[g] for g in cohort_list]
        total = sum(p)
        centered = [w1[i] * (G[i] == g) - share[g] for i in 1:n, g in cohort_list]
        att = first.(by_cohort)
        (sum(att .* p) / total, reduce(hcat, last.(by_cohort)) * (p ./ total) .+
            (centered ./ total .- sum(centered; dims = 2) * (p' ./ total^2)) * att)
    end
    over_times = averaged(by_time)
    overall = (; aggregation = ["simple", "event time", "cohort", "calendar"],
        summary([simple, over_events, over_cohorts, over_times])...)

    sample = (; units = n, clusters = C, times = time_list, cohorts = groups,
        never_treated,
        control = control == "never" ? "never treated" : control == "notyet" ? "not yet treated" :
            "never treated and treated last ($latest), before their treatment",
        base_period, method, anticipation, dropped = length(always), unbalanced = gaps,
        panel = repeated ? "unbalanced" : "balanced", draws = bootstrap ? draws : 0, uniform,
        critical = (; cells = cell_critical, effects = event_critical, cohorts = cohort_critical,
            calendar = calendar_critical, overall = critical))
    CallawaySantAnnaResult(sample, cells, effects, cohort_table, calendar, overall, pretrends)
end

"""
The 2x2 difference in differences of the changes `dy` of the treated (`D` = 1) and the control units, with the
covariates `X` (a constant first) and the weights `w`, by `method`, and its influence function (as DRDID's
drdid_panel, std_ipw_did_panel and reg_did_panel).
"""
function two_by_two(dy, D, X, w, method)
    n = length(D)
    w = w ./ mean(w)
    out = zeros(n)
    if method != "ipw"
        c = D .== 0
        beta = (X[c, :]' * (w[c] .* X[c, :])) \ (X[c, :]' * (w[c] .* dy[c]))
        out = X * beta
        XpX = X' * ((w .* (1 .- D)) .* X) / n
        asy_ols = ((w .* (1 .- D) .* (dy .- out)) .* X) / XpX
    end
    if method == "reg"
        w_treat = w .* D
        eta_treat = mean(w_treat .* dy) / mean(w_treat)
        eta_cont = mean(w_treat .* out) / mean(w_treat)
        inf_treat = (w_treat .* dy .- w_treat .* eta_treat) ./ mean(w_treat)
        inf_control = (w_treat .* out .- w_treat .* eta_cont .+ asy_ols * (vec(w_treat' * X) / n)) ./ mean(w_treat)
        return (eta_treat - eta_cont, inf_treat .- inf_control)
    end
    ps = min.(logit_fit(X, D, w), 1 - 1e-6)
    trim = [D[i] == 1 || ps[i] < 0.995 for i in 1:n]
    w_treat = trim .* w .* D
    w_cont = trim .* w .* ps .* (1 .- D) ./ (1 .- ps)
    residual = dy .- out
    eta_treat = mean(w_treat .* residual) / mean(w_treat)
    eta_cont = mean(w_cont .* residual) / mean(w_cont)
    asy_ps = ((w .* (D .- ps)) .* X) * (inv(cholesky(Symmetric(X' * ((ps .* (1 .- ps) .* w) .* X)))) * n)
    inf_treat = w_treat .* residual .- w_treat .* eta_treat
    inf_cont = w_cont .* residual .- w_cont .* eta_cont .+ asy_ps * (vec((w_cont .* (residual .- eta_cont))' * X) / n)
    if method == "dr"
        inf_treat .-= asy_ols * (vec(w_treat' * X) / n)
        inf_cont .-= asy_ols * (vec(w_cont' * X) / n)
    end
    (eta_treat - eta_cont, inf_treat ./ mean(w_treat) .- inf_cont ./ mean(w_cont))
end

"""
The 2x2 difference in differences of repeated cross sections: the outcomes `y` before (`post` = 0) and after
(`post` = 1) of the treated (`D` = 1) and the control units, with the covariates `X` (a constant first) and the weights
`w`, by `method`, and its influence function, an observation a row (as DRDID's drdid_rc, std_ipw_did_rc and
reg_did_rc).
"""
function two_by_two_rc(y, post, D, X, w, method)
    n = length(D)
    w = w ./ mean(w)
    # The weighted mean of r under the weights v, its denominator and its influence
    function weighed(v, r)
        m = mean(v)
        a = mean(v .* r) / m
        (a, m, v .* (r .- a) ./ m)
    end
    # The weighted least squares fit of y on X over the rows `which`, everywhere, and its asymptotic linear representation
    function ols(which)
        v = w .* which
        beta = (X' * (v .* X)) \ (X' * (v .* y))
        fit = X * beta
        (fit, ((v .* (y .- fit)) .* X) / (X' * (v .* X) / n))
    end
    if method != "ipw"
        control_pre, asy_pre = ols((D .== 0) .& (post .== 0))
        control_post, asy_post = ols((D .== 0) .& (post .== 1))
    end
    if method == "reg"
        a_pre, _, inf_pre = weighed(w .* D .* (1 .- post), y)
        a_post, _, inf_post = weighed(w .* D .* post, y)
        w_cont = w .* D
        change = control_post .- control_pre
        a_cont, m_cont, _ = weighed(w_cont, change)
        M1 = vec(w_cont' * X) / n
        inf_control = (w_cont .* (change .- a_cont) .+ asy_post * M1 .- asy_pre * M1) ./ m_cont
        return ((a_post - a_pre) - a_cont, (inf_post .- inf_pre) .- inf_control)
    end
    ps = min.(logit_fit(X, D, w), 1 - 1e-6)
    trim = [D[i] == 1 || ps[i] < 0.995 for i in 1:n]
    w_treat_pre = trim .* w .* D .* (1 .- post)
    w_treat_post = trim .* w .* D .* post
    odds = trim .* w .* ps .* (1 .- D) ./ (1 .- ps)
    w_cont_pre = odds .* (1 .- post)
    w_cont_post = odds .* post
    asy_ps = ((w .* (D .- ps)) .* X) * (inv(cholesky(Symmetric(X' * ((ps .* (1 .- ps) .* w) .* X)))) * n)
    # Inverse probability weighting on the outcomes, doubly robust on what the control units' regressions leave
    residual = method == "dr" ? y .- (post .* control_post .+ (1 .- post) .* control_pre) : y
    a_treat_pre, m_treat_pre, inf_treat_pre = weighed(w_treat_pre, residual)
    a_treat_post, m_treat_post, inf_treat_post = weighed(w_treat_post, residual)
    a_cont_pre, m_cont_pre, inf_cont_pre = weighed(w_cont_pre, residual)
    a_cont_post, m_cont_post, inf_cont_post = weighed(w_cont_post, residual)
    M2_pre = vec((w_cont_pre .* (residual .- a_cont_pre))' * X) / n / m_cont_pre
    M2_post = vec((w_cont_post .* (residual .- a_cont_post))' * X) / n / m_cont_post
    inf_treat = inf_treat_post .- inf_treat_pre
    inf_cont = inf_cont_post .- inf_cont_pre .+ asy_ps * (M2_post .- M2_pre)
    att = (a_treat_post - a_treat_pre) - (a_cont_post - a_cont_pre)
    method == "ipw" && return (att, inf_treat .- inf_cont)

    # The doubly robust estimator's terms from the treated units' regressions
    treated_pre, asy_treated_pre = ols((D .== 1) .& (post .== 0))
    treated_post, asy_treated_post = ols((D .== 1) .& (post .== 1))
    inf_treat .+= asy_post * (-vec(w_treat_post' * X) / n / m_treat_post) .+ asy_pre * (-vec(w_treat_pre' * X) / n / m_treat_pre)
    inf_cont .+= asy_post * (-vec(w_cont_post' * X) / n / m_cont_post) .+ asy_pre * (-vec(w_cont_pre' * X) / n / m_cont_pre)
    w_d = trim .* w .* D
    w_d1 = w_d .* post
    w_d0 = w_d .* (1 .- post)
    gap_post = treated_post .- control_post
    gap_pre = treated_pre .- control_pre
    a_d_post, m_d, inf_d_post = weighed(w_d, gap_post)
    a_d1, m_d1, inf_d1 = weighed(w_d1, gap_post)
    a_d_pre, _, inf_d_pre = weighed(w_d, gap_pre)
    a_d0, m_d0, inf_d0 = weighed(w_d0, gap_pre)
    att += (a_d_post - a_d1) - (a_d_pre - a_d0)
    inf_effect = (inf_d_post .- inf_d1) .- (inf_d_pre .- inf_d0)
    moment_post = vec((w_d ./ m_d .- w_d1 ./ m_d1)' * X) / n
    moment_pre = vec((w_d ./ m_d .- w_d0 ./ m_d0)' * X) / n
    inf_or = (asy_treated_post .- asy_post) * moment_post .- (asy_treated_pre .- asy_pre) * moment_pre
    (att, inf_treat .- inf_cont .+ inf_effect .+ inf_or)
end

"""The fitted probabilities of a weighted logit of `D` on `X`, by Newton's method."""
function logit_fit(X, D, w)
    beta = zeros(size(X, 2))
    p = fill(0.5, length(D))
    for _ in 1:100
        p = 1 ./ (1 .+ exp.(-(X * beta)))
        step = cholesky(Symmetric(X' * ((w .* p .* (1 .- p)) .* X))) \ (X' * (w .* (D .- p)))
        beta .+= step
        maximum(abs, step) < 1e-12 && break
    end
    1 ./ (1 .+ exp.(-(X * beta)))
end

function report(io::IO, result::CallawaySantAnnaResult)
    s = result.sample
    println(io, "\nCallaway and Sant'Anna: ", s.units, " units in ", s.clusters, " clusters, times ", first(s.times), " to ",
        last(s.times), "; cohorts ", join(s.cohorts, ", "), "; control: ", s.control, "; base period ", s.base_period,
        "; method ", s.method, s.anticipation > 0 ? "; anticipation $(s.anticipation)" : "",
        s.dropped > 0 ? "; $(s.dropped) units treated from the first time left out" : "",
        s.unbalanced == 0 ? "" : s.panel == "unbalanced" ?
            "; $(s.unbalanced) units not seen at every time kept, the cells as repeated cross sections" :
            "; $(s.unbalanced) units not seen at every time left out")
    c = s.critical
    s.draws == 0 ? println(io, "Standard errors from the influence functions, clustered; 95% bands pointwise") :
        println(io, "Standard errors from the multiplier bootstrap, ", s.draws, " draws by cluster; 95% bands ",
            s.uniform ? "simultaneous, critical values $(format(c.cells)) (cells), $(format(c.effects)) (event times), " *
            "$(format(c.cohorts)) (cohorts), $(format(c.calendar)) (times), pointwise $(format(c.overall)) (overall)" : "pointwise")
    for (title, key, table) in (("event", :event_time, result.effects), ("cohort", :cohort, result.cohorts),
            ("time", :time, result.calendar), ("overall", :aggregation, result.overall))
        @printf(io, "%10s   %9s %8s %8s %7s %9s %9s\n", title, "effect", "se", "t", "p", "low", "high")
        for j in eachindex(getproperty(table, key))
            @printf(io, "%10s   %9s %8s %8s %7s %9s %9s\n", getproperty(table, key)[j], format(table.effect[j]),
                format(table.se[j]), format(table.t[j]), format(table.p[j]), format(table.low[j]), format(table.high[j]))
        end
    end
    p = result.pretrends
    degrees = p.test == "F" ? "F($(p.df1), $(p.df2))" : "chi2($(p.df1))"
    ismissing(p.statistic) ? println(io, "Pre-trends: not tested (", p.df1 == 0 ? "no cells before treatment" :
        "the covariance of the $(p.df1) cells before treatment is singular, with $(s.clusters) clusters", ")") :
        println(io, "Pre-trends: ", degrees, " = ", format(p.statistic), ", p = ", format(p.p))
end

"""
    write_results(result::CallawaySantAnnaResult, folder)

Writes the cells (`cs-cells.csv`), the effects by event time (`cs-effects.csv`), by cohort (`cs-cohorts.csv`), by
calendar time (`cs-calendar.csv`), the overall averages (`cs-overall.csv`) and the test of pre-trends
(`cs-pretrends.csv`) to `folder`; returns the files.
"""
function write_results(result::CallawaySantAnnaResult, folder::AbstractString)
    mkpath(folder)
    [CSV.write(joinpath(folder, name), table) for (name, table) in result_tables(result)]
end

result_tables(result::CallawaySantAnnaResult) = ["cs-cells.csv" => result.cells, "cs-effects.csv" => result.effects,
    "cs-cohorts.csv" => result.cohorts, "cs-calendar.csv" => result.calendar, "cs-overall.csv" => result.overall,
    "cs-pretrends.csv" => map(x -> [x], result.pretrends)]
