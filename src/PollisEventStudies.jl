"""
    PollisEventStudies

Intraday event studies of SEC filings and of the FOMC's releases, on the folders Pollis keeps for them.

A study of filings has `study.toml` (the benchmark, the model of normal returns and the periods), `events.csv` (the
filings) and `prices/` (a file per filing with the minute bars of the stock and of the benchmark). For each period: the
average and cumulative average abnormal returns a minute at a time with their cross-sectional t, the Boehmer, Musumeci
and Poulsen (1991) test on standardized CARs, and abnormal volume.

A study of FOMC releases (`source = "fomc"` in `study.toml`) has `events.csv` (the releases), `instruments.csv` (each
instrument's contract, unit and scale at each release) and `prices/` (a file per release with the bars of each
contract). For each period: each instrument's move at each release, the mean move with its cross-sectional t, and the
mean move a minute at a time.

A panel event study (`source = "panel"`) has a table of units over time, first treated at different times: the effects
by event time from two-way fixed effects and from Sun and Abraham (2021), each cohort's, and the tests of pre-trends.
`panel_event_study` takes any Tables.jl table, as does `callaway_santanna`, for Callaway and Sant'Anna's (2021) effects
of each cohort at each time and their aggregations.

    results = run_study("~/.pollis/event-studies/my-study")
"""
module PollisEventStudies

using CSV, Tables, TOML, Dates, Statistics, Printf, LinearAlgebra
import Random
using Distributions: TDist, FDist, Normal, Chisq, ccdf, quantile

export read_study, analyse, report, write_results, run_study, instrument_move, panel_event_study, callaway_santanna

# A minute's bar is complete a little after it ends
const SETTLE = Minute(2)
# The fewest minute returns alpha, beta and the residual standard deviation are estimated on
const MIN_RETURNS = 30
# The fewest regular-session minutes before the filing abnormal volume is measured against
const MIN_VOLUME_MINUTES = 10

"""Minutes before and after the filing; `next_open` also follows a filing outside the regular session from the next open."""
struct Period
    name::String
    before::Int
    after::Int
    next_open::Bool
end

"""
A filing: its form and 8-K items, when EDGAR accepted it and the event's time (UTC), which is the earlier of the
newswire release and the acceptance when the study takes the first public time; whether the event came during the
regular session, and the next open after it.
"""
struct Event
    accession::String
    cik::String
    name::String
    ticker::String
    form::String
    items::Vector{String}
    accepted::DateTime
    time::DateTime
    during_session::Union{Bool, Nothing}
    next_open::Union{DateTime, Nothing}
    complete::Bool
end

struct Study
    folder::String
    title::String
    benchmark::String
    market_model::Bool
    estimation_days::Int
    periods::Vector{Period}
    events::Vector{Event}
    # What the events are grouped by: "form", "item", "session", "ticker", or "" for no groups
    group_by::String
end

"""A symbol's minute bars, by when each minute ends."""
struct Series
    ends::Vector{DateTime}
    close::Vector{Float64}
    volume::Vector{Float64}
    regular::Vector{Bool}
end

# The times in the files are UTC, in ISO 8601 with a Z
utc(text) = DateTime(replace(text, "Z" => ""))

"""
    read_study(folder = pwd())

The study in `folder`: its settings and the events with a ticker and an acceptance time.
"""
function read_study(folder::AbstractString = pwd())
    settings = TOML.parsefile(joinpath(folder, "study.toml"))
    get(settings, "source", "filings") == "panel" && return read_panel_study(folder, settings)
    fomc = get(settings, "source", "filings") == "fomc"
    periods = [Period(get(p, "name", "Period $k"), get(p, "before", 30), get(p, "after", 120), get(p, "next_open", false))
               for (k, p) in enumerate(get(settings, "periods", Dict{String, Any}[]))]
    isempty(periods) && push!(periods, fomc ? Period("Tight", 10, 20, false) : Period("Default", 30, 120, true))
    fomc && return read_fomc_study(folder, settings, periods)
    events = Event[]
    path = joinpath(folder, "events.csv")
    if isfile(path)
        file = CSV.File(path; types = String)
        # An events.csv written before Pollis looked for newswire releases has no event_time: the acceptance is the event's time
        has_time = :event_time in propertynames(file)
        has_items = :items in propertynames(file)
        for row in Tables.rows(file)
            (ismissing(row.ticker) || ismissing(row.accepted)) && continue
            during = coalesce(row.during_session, "")
            time = has_time && !ismissing(row.event_time) ? row.event_time : row.accepted
            # The 8-K items, separated by spaces
            items = has_items && !ismissing(row.items) ? String.(split(row.items)) : String[]
            push!(events, Event(row.accession, coalesce(row.cik, ""), coalesce(row.name, ""), row.ticker, coalesce(row.form, ""),
                items, utc(row.accepted), utc(time), during == "true" ? true : during == "false" ? false : nothing,
                ismissing(row.next_open) ? nothing : utc(row.next_open), coalesce(row.status, "") == "complete"))
        end
    end
    Study(folder, get(settings, "title", basename(folder)), get(settings, "benchmark", "^GSPC"),
        get(settings, "model", "marketAdjusted") == "marketModel", get(settings, "estimation_days", 5), periods, events,
        get(settings, "group_by", ""))
end

"""The minute bars of an event's prices file, by symbol."""
read_prices(study::Study, event::Event) = read_bars(joinpath(study.folder, "prices", event.accession * ".csv"), 1)

"""The bars of a prices file, by symbol, each spanning `minutes` from its time."""
function read_bars(path, minutes)
    series = Dict{String, Series}()
    isfile(path) || return series
    rows = Tables.rowtable(CSV.File(path; types = Dict(:time => String, :session => String, :symbol => String)))
    for row in sort(rows; by = row -> row.time)
        bars = get!(() -> Series(DateTime[], Float64[], Float64[], Bool[]), series, row.symbol)
        push!(bars.ends, utc(row.time) + Minute(minutes))
        push!(bars.close, row.close)
        push!(bars.volume, coalesce(row.volume, 0))
        push!(bars.regular, row.session == "regular")
    end
    series
end

"""The last price at or before time t; nothing before the first bar."""
function price_at(bars::Series, t::DateTime)
    k = searchsortedlast(bars.ends, t)
    k == 0 ? nothing : bars.close[k]
end

"""The bar that ends at time t; nothing when there is none (a minute without trades has no bar)."""
function bar_at(bars::Series, t::DateTime)
    k = searchsortedlast(bars.ends, t)
    k > 0 && bars.ends[k] == t ? k : nothing
end

"""How many of the bars' minutes end after from and at or before to."""
count_between(bars::Series, from::DateTime, to::DateTime) = max(0, searchsortedlast(bars.ends, to) - searchsortedlast(bars.ends, from))

"""
How many minutes ending after from and at or before to the stock or the benchmark has a bar in: every minute of the
regular session, whether the stock traded in it or not, and outside it those the stock traded in.
"""
function window_minutes(stock::Series, benchmark::Series, from::DateTime, to::DateTime)
    span(bars) = bars.ends[searchsortedlast(bars.ends, from) + 1:searchsortedlast(bars.ends, to)]
    length(union(span(stock), span(benchmark)))
end

"""
The stock's log returns from one trade to the next in the regular sessions of the `days` trading days before `start`,
the benchmark's over the same minutes, how many minutes each spans and when each ends. A stock that does not trade every minute moves
over the minutes it did not trade in at its next trade, which is measured against the benchmark's move over all of
them (thin trading, as `cum_periods` in eventstudy2 and in EventStudy.jl), rather than a return of zero in each.
"""
function estimation_returns(stock::Series, benchmark::Series, start::DateTime, days)
    # The benchmark trades only in the regular session: its minutes are the session's, all in one UTC day
    n = searchsortedlast(benchmark.ends, start - Millisecond(1))
    ends = benchmark.ends[1:n]
    chosen = Set(unique(Date.(ends))[max(1, end - days + 1):end])
    stock_returns, benchmark_returns, minutes, times = Float64[], Float64[], Int[], DateTime[]
    # The benchmark's minute of the stock's last trade in the session, and the stock's bar then
    last_minute, last_bar = 0, 0
    for k in 1:n
        Date(ends[k]) in chosen || continue
        # A gap in the benchmark's minutes starts another session: no return spans it
        (k == 1 || ends[k] - ends[k - 1] != Minute(1)) && (last_minute = 0)
        bar = bar_at(stock, ends[k])
        bar === nothing && continue
        if last_minute > 0
            push!(stock_returns, log(stock.close[bar] / stock.close[last_bar]))
            push!(benchmark_returns, log(benchmark.close[k] / benchmark.close[last_minute]))
            push!(minutes, k - last_minute)
            push!(times, ends[k])
        end
        last_minute, last_bar = k, bar
    end
    stock_returns, benchmark_returns, minutes, times
end

"""
Alpha (per minute of the regular session), beta and the residual standard deviation of a minute; nothing with too
few returns. A return over g minutes has g times a minute's alpha and g times its variance: the returns are weighed
by 1 / g (weighted least squares), which is ordinary least squares when the stock trades every minute.
"""
function fit(stock_returns, benchmark_returns, minutes)
    length(stock_returns) < MIN_RETURNS && return nothing
    # The normal equations of r / sqrt(g) = alpha sqrt(g) + beta m / sqrt(g)
    g, m, r = minutes, benchmark_returns, stock_returns
    gg, gm, mm, gr, mr = sum(g), sum(m), sum(m .^ 2 ./ g), sum(r), sum(m .* r ./ g)
    determinant = gg * mm - gm^2
    determinant <= 0 && return nothing
    alpha, beta = (mm * gr - gm * mr) / determinant, (gg * mr - gm * gr) / determinant
    residuals = (r .- alpha .* g .- beta .* m) ./ sqrt.(g)
    (alpha = alpha, beta = beta, sigma = sqrt(sum(residuals .^ 2) / (length(residuals) - 2)), residuals = residuals)
end

"""
An event's cumulative abnormal returns over a period, in per cent from the start of the period, a value a minute from
the filing and, for a filing outside the regular session in a period that follows the next open, from the next open;
and at the end of the period (for such a filing, the period's minutes after the next open) its CAR, standardized CAR,
abnormal volume and beta. Nothing without the prices at its start.
"""
function analyse_event(study::Study, period::Period, event::Event, now::DateTime)
    prices = read_prices(study, event)
    stock, benchmark = get(prices, event.ticker, nothing), get(prices, study.benchmark, nothing)
    (stock === nothing || benchmark === nothing) && return nothing
    start = event.time - Minute(period.before)
    stock_start, benchmark_start = price_at(stock, start), price_at(benchmark, start)
    (stock_start === nothing || benchmark_start === nothing) && return nothing

    stock_returns, benchmark_returns, minutes, times = estimation_returns(stock, benchmark, start, study.estimation_days)
    model = study.market_model ? fit(stock_returns, benchmark_returns, minutes) : nothing
    residuals = model !== nothing ? model.residuals :
        length(stock_returns) >= MIN_RETURNS ? (stock_returns .- benchmark_returns) ./ sqrt.(minutes) : nothing
    sigma = model !== nothing ? model.sigma : residuals !== nothing ? std(residuals) : nothing
    known = event.complete ? typemax(DateTime) : now - SETTLE

    function car(t)
        price, level = price_at(stock, t), price_at(benchmark, t)
        (t > known || price === nothing || level === nothing) && return missing
        market = log(level / benchmark_start)
        normal = model === nothing ? market : model.alpha * count_between(benchmark, start, t) + model.beta * market
        100 * (log(price / stock_start) - normal)
    end

    from_filing = Union{Missing, Float64}[car(event.time + Minute(m)) for m in -period.before:period.after]
    follows_open = period.next_open && event.during_session === false && event.next_open !== nothing
    from_open = follows_open ? Union{Missing, Float64}[car(event.next_open + Minute(m)) for m in 0:period.after] : nothing

    # A filing outside the regular session followed from the next open ends the period after the next open, where the
    # market first trades on it: its CAR from the filing would end while the market is closed
    finish = (follows_open ? event.next_open : event.time) + Minute(period.after)
    final = follows_open ? from_open[end] : from_filing[end]
    scar = sigma === nothing || ismissing(final) ? missing :
        final / 100 / (sigma * sqrt(max(1, window_minutes(stock, benchmark, start, finish))))

    # The mean log(1 + volume) of the regular-session minutes after the event, or after the next open, less that of
    # the regular-session minutes before the event
    from = follows_open ? event.next_open : event.time
    window = [log1p(stock.volume[k]) for k in eachindex(stock.ends) if stock.regular[k] && from < stock.ends[k] <= min(from + Minute(period.after), known)]
    baseline = [log1p(stock.volume[k]) for k in eachindex(stock.ends) if stock.regular[k] && stock.ends[k] <= event.time]
    volume = isempty(window) || length(baseline) < MIN_VOLUME_MINUTES ? missing : mean(window) - mean(baseline)

    (accession = event.accession, ticker = event.ticker, name = event.name, from_filing = from_filing, from_open = from_open,
        car = final, scar = scar, abnormal_volume = volume, beta = model === nothing ? missing : model.beta,
        # The standardized abnormal returns of the estimation days, by when each ends, and the period, for the tests
        # of Kolari and Pynnonen
        sar_times = times, sars = sigma === nothing ? Float64[] : residuals ./ sigma, window = (start, finish))
end

"""The cross-sectional t of a sample: its mean over its standard error; missing with fewer than two."""
tstat(values) = length(values) < 2 ? missing : mean(values) / (std(values) / sqrt(length(values)))

"""The two-sided p-value of a t with `df` degrees of freedom; missing without one."""
pvalue(t, df) = ismissing(t) || !isfinite(t) || df < 1 ? missing : 2 * ccdf(TDist(df), abs(t))

"""
The mean cross-correlation of the events' standardized abnormal returns (Kolari and Pynnonen 2010), over every pair
of events: over the minutes of their estimation days in common for a pair whose periods overlap, and zero for a pair
whose periods do not, whose abnormal returns come from different minutes. Missing with fewer than two events.
"""
function cross_correlation(results)
    n = length(results)
    n < 2 && return missing
    total = 0.0
    for i in 1:n - 1, j in i + 1:n
        a, b = results[i], results[j]
        (a.window[2] <= b.window[1] || b.window[2] <= a.window[1]) && continue
        at = Dict(zip(a.sar_times, a.sars))
        pairs = [(at[time], b.sars[k]) for (k, time) in enumerate(b.sar_times) if haskey(at, time)]
        length(pairs) < MIN_RETURNS && continue
        r = cor(first.(pairs), last.(pairs))
        isfinite(r) && (total += r)
    end
    total / (n * (n - 1) / 2)
end

"""
How much the t of the standardized CARs shrinks for their cross-correlation r (Kolari and Pynnonen 2010):
sqrt((1 - r) / (1 + (n - 1) r)); 1 when the events' periods do not overlap.
"""
function correlation_factor(r, n)
    (ismissing(r) || 1 + (n - 1) * r <= 0) && return missing
    sqrt((1 - r) / (1 + (n - 1) * r))
end

"""
The generalized rank test of Kolari and Pynnonen (2011) on the events' CARs, robust to returns far from normal, to
the volatility an event brings and to the events' cross-correlation: each event's standardized CAR, restandardized by
their cross-sectional standard deviation, is ranked among the standardized abnormal returns of its estimation days.
Its t and degrees of freedom; missing with fewer than two events.
"""
function generalized_rank(results)
    n = length(results)
    missing_test = (t = missing, df = 0)
    n < 2 && return missing_test
    scale = std([r.scar for r in results])
    (isfinite(scale) && scale > 0) || return missing_test
    # The mean standardized rank at each time, counted back from the event (1) through the estimation days, over the
    # events with an abnormal return then
    sums, counts = Float64[], Int[]
    for result in results
        values = [result.sars; result.scar / scale]
        ranks = invperm(sortperm(values))
        for (k, rank) in enumerate(reverse(ranks))
            k > length(sums) && (push!(sums, 0.0); push!(counts, 0))
            sums[k] += rank / (length(values) + 1) - 0.5
            counts[k] += 1
        end
    end
    times = length(sums)
    times < 3 && return missing_test
    means = sums ./ counts
    spread = sqrt(sum(counts ./ n .* means .^ 2) / times)
    z = means[1] / spread
    times - 1 - z^2 <= 0 && return missing_test
    (t = z * sqrt((times - 2) / (times - 1 - z^2)), df = times - 2)
end

"""The average of the events' curves a minute at a time, over the events known at that minute, and its t."""
function average(curves, minutes)
    n = length(minutes)
    counts = zeros(Int, n)
    means = Vector{Union{Missing, Float64}}(missing, n)
    ts = Vector{Union{Missing, Float64}}(missing, n)
    for k in 1:n
        values = Float64[curve[k] for curve in curves if !ismissing(curve[k])]
        counts[k] = length(values)
        isempty(values) && continue
        means[k] = mean(values)
        ts[k] = tstat(values)
    end
    (minute = collect(minutes), events = counts, mean = means, t = ts)
end

"""Each minute's abnormal return: the change in the cumulative one."""
increments(curve) = Union{Missing, Float64}[k == 1 ? missing : curve[k] - curve[k - 1] for k in eachindex(curve)]

"""The AAR and CAAR of the curves a minute at a time, with their t, as a table."""
function averages(curves, minutes)
    caar = average(curves, minutes)
    aar = average(increments.(curves), minutes)
    (minute = caar.minute, events = caar.events, aar = aar.mean, aar_t = aar.t, caar = caar.mean, caar_t = caar.t)
end

known(values) = Float64[value for value in values if !ismissing(value)]

"""
The events' results over a period: each event's (`events`, a table), each event's CAR a minute at a time from the
filing (`curves`), the AAR and CAAR from the filing (`from_filing`) and from the next open (`from_open`), and the tests
at the end of the period (`tests`).
"""
function summarise(period::Period, results)
    opens = [result.from_open for result in results if result.from_open !== nothing]
    cars, scars, volumes = known(r.car for r in results), known(r.scar for r in results), known(r.abnormal_volume for r in results)
    mean_or_missing(values) = isempty(values) ? missing : mean(values)
    standardized = [r for r in results if !ismissing(r.scar)]
    correlation = cross_correlation(standardized)
    factor = correlation_factor(correlation, length(standardized))
    bmp_t = tstat(scars)
    bmp_kp_t = ismissing(bmp_t) || ismissing(factor) ? missing : bmp_t * factor
    grank = generalized_rank(standardized)
    (events = (accession = [r.accession for r in results], ticker = [r.ticker for r in results], name = [r.name for r in results],
            car = [r.car for r in results], scar = [r.scar for r in results], abnormal_volume = [r.abnormal_volume for r in results],
            beta = [r.beta for r in results]),
        curves = [r.from_filing for r in results],
        from_filing = isempty(results) ? nothing : averages([r.from_filing for r in results], -period.before:period.after),
        from_open = isempty(opens) ? nothing : averages(opens, 0:period.after),
        open_events = length(opens),
        tests = (car = mean_or_missing(cars), car_t = tstat(cars), car_p = pvalue(tstat(cars), length(cars) - 1), car_events = length(cars),
            bmp = mean_or_missing(scars), bmp_t = bmp_t, bmp_p = pvalue(bmp_t, length(scars) - 1), bmp_events = length(scars),
            correlation = correlation, bmp_kp_t = bmp_kp_t, bmp_kp_p = pvalue(bmp_kp_t, length(scars) - 1),
            grank_t = grank.t, grank_p = pvalue(grank.t, grank.df),
            abnormal_volume = mean_or_missing(volumes), abnormal_volume_t = tstat(volumes),
            abnormal_volume_p = pvalue(tstat(volumes), length(volumes) - 1), abnormal_volume_events = length(volumes)))
end

"""
The groups an event is in: its form, each of its 8-K items (a filing with several is in the group of each), whether
it came during the regular session, its ticker, or the name or names a function of the event gives.
"""
function groups_of(event::Event, group_by)
    group_by isa Function && return string.(vcat(group_by(event)))
    group_by == "form" && return [event.form]
    group_by == "item" && return isempty(event.items) ? ["No item"] : ["Item " * item for item in event.items]
    group_by == "session" && return [event.during_session === nothing ? "Session unknown" :
        event.during_session ? "During the session" : "Outside the session"]
    group_by == "ticker" && return [event.ticker]
    throw(ArgumentError("group_by is \"form\", \"item\", \"session\", \"ticker\" or a function of the event, not \"$group_by\""))
end

"""What the report calls the groups."""
grouping(group_by) = group_by isa Function ? "group" : get(Dict("form" => "form", "item" => "8-K item",
    "session" => "session", "ticker" => "ticker"), group_by, group_by)

"""
    analyse(study; group_by = study.group_by, now = Dates.now(UTC))

For each period of the study: its events' results (`events`, a table), each event's CAR a minute at a time from
the filing (`curves`), the AAR and CAAR from the filing
(`from_filing`, a table a minute a row) and from the next open (`from_open`, or nothing), and the tests at the end of
the period (`tests`). Events still capturing prices count up to the last complete minute before `now`.

With `group_by` ("form", "item", "session", "ticker", or a function of the event that gives a group's name or
several), the same for each group of events (`groups`, each with its `name`), as EventStudy.jl's groups do.
"""
function analyse(study::Study; group_by = study.group_by, now::DateTime = Dates.now(UTC))
    grouped = !(group_by isa AbstractString && isempty(group_by))
    map(study.periods) do period
        analysed = [(event, result) for (event, result) in ((event, analyse_event(study, period, event, now)) for event in study.events) if result !== nothing]
        results = last.(analysed)
        names = grouped ? sort(unique(reduce(vcat, (groups_of(event, group_by) for (event, _) in analysed); init = String[]))) : String[]
        groups = [merge((name = name,), summarise(period, [result for (event, result) in analysed if name in groups_of(event, group_by)]))
                  for name in names]
        merge((period = period,), summarise(period, results), (grouped_by = grouped ? grouping(group_by) : "", groups = groups))
    end
end

format(x) = ismissing(x) || !isfinite(x) ? "" : @sprintf("%.3f", x)

function report_table(io::IO, table, marks)
    @printf(io, "%8s %7s %10s %8s %10s %8s\n", "minute", "events", "AAR %", "t", "CAAR %", "t")
    for k in eachindex(table.minute)
        table.minute[k] in marks || continue
        @printf(io, "%8d %7d %10s %8s %10s %8s\n", table.minute[k], table.events[k], format(table.aar[k]), format(table.aar_t[k]), format(table.caar[k]), format(table.caar_t[k]))
    end
end

"""
    report([io,] results)

Prints, for each period, the AAR and CAAR at some minutes from the filing and from the next open, and the tests at the
end of the period.
"""
report(results) = report(stdout, results)
function report(io::IO, results)
    for result in results
        period = result.period
        println(io, "\n", period.name, ": ", period.before, " minutes before the filing to ", period.after, " after; ",
            length(result.events.accession), " events with prices")
        result.from_filing === nothing && continue
        marks = filter(m -> -period.before <= m <= period.after, [-period.before, -15, -5, 0, 1, 5, 15, 30, 60, 120, 240, period.after])
        println(io, "\nMinutes from the filing")
        report_table(io, result.from_filing, marks)
        if result.from_open !== nothing
            println(io, "\nMinutes from the next open, ", result.open_events, " events")
            report_table(io, result.from_open, marks)
        end
        tests = result.tests
        @printf(io, "\nAt +%d minutes: CAR %s%% (t = %s, p = %s, %d events)\n", period.after, format(tests.car), format(tests.car_t), format(tests.car_p), tests.car_events)
        if tests.bmp_events == 0
            println(io, "Boehmer, Musumeci and Poulsen: needs the estimation days, which only the market model captures")
        else
            @printf(io, "Boehmer, Musumeci and Poulsen: mean standardized CAR %s (t = %s, p = %s, %d events)\n", format(tests.bmp), format(tests.bmp_t), format(tests.bmp_p), tests.bmp_events)
            @printf(io, "  adjusted for cross-correlation (Kolari and Pynnonen 2010, mean correlation %s): t = %s, p = %s\n",
                format(tests.correlation), format(tests.bmp_kp_t), format(tests.bmp_kp_p))
            @printf(io, "Generalized rank test (Kolari and Pynnonen 2011): t = %s, p = %s\n", format(tests.grank_t), format(tests.grank_p))
        end
        @printf(io, "Abnormal volume: %s (t = %s, p = %s, %d events)\n", format(tests.abnormal_volume), format(tests.abnormal_volume_t), format(tests.abnormal_volume_p), tests.abnormal_volume_events)
        isempty(result.groups) && continue
        println(io, "\nBy ", result.grouped_by, ", at +", period.after, " minutes")
        @printf(io, "%-24s %7s %10s %8s %8s %8s %8s %10s %8s\n", result.grouped_by, "events", "CAR %", "t", "BMP t", "KP t", "GRANK t", "volume", "t")
        for group in groups_table(result)
            @printf(io, "%-24s %7d %10s %8s %8s %8s %8s %10s %8s\n", first(group.group, 24), group.events, format(group.car), format(group.car_t),
                format(group.bmp_t), format(group.bmp_kp_t), format(group.grank_t), format(group.abnormal_volume), format(group.abnormal_volume_t))
        end
    end
end

"""A period's groups a row, with their tests at the end of the period."""
groups_table(result) = [(group = group.name, events = length(group.events.accession), car = group.tests.car, car_t = group.tests.car_t,
    car_p = group.tests.car_p, bmp = group.tests.bmp, bmp_t = group.tests.bmp_t, bmp_p = group.tests.bmp_p,
    correlation = group.tests.correlation, bmp_kp_t = group.tests.bmp_kp_t, bmp_kp_p = group.tests.bmp_kp_p,
    grank_t = group.tests.grank_t, grank_p = group.tests.grank_p, abnormal_volume = group.tests.abnormal_volume,
    abnormal_volume_t = group.tests.abnormal_volume_t, abnormal_volume_p = group.tests.abnormal_volume_p) for group in result.groups]

"""The tables of a period's groups one after the other, each row with its group's name."""
stacked(result, table) = [merge((group = group.name,), row) for group in result.groups if getproperty(group, table) !== nothing
                          for row in Tables.rowtable(getproperty(group, table))]

slug(name) = strip(replace(lowercase(name), r"[^a-z0-9]+" => "-"), '-')

"""
    write_results(results, folder)

Writes, for each period, its events (`<period>-events.csv`) and its AAR and CAAR a minute a row
(`<period>-from-the-filing.csv`, `<period>-from-the-next-open.csv`) to `folder` and, when the events are grouped, the
groups' tests (`<period>-groups.csv`, a group a row) and their AAR and CAAR (`<period>-groups-from-the-filing.csv`,
`<period>-groups-from-the-next-open.csv`, with each row's group); returns the files.
"""
function write_results(results, folder::AbstractString)
    mkpath(folder)
    files = String[]
    save(name, table) = push!(files, CSV.write(joinpath(folder, name), table))
    for result in results
        name = slug(result.period.name)
        save(name * "-events.csv", result.events)
        result.from_filing === nothing || save(name * "-from-the-filing.csv", result.from_filing)
        result.from_open === nothing || save(name * "-from-the-next-open.csv", result.from_open)
        isempty(result.groups) && continue
        save(name * "-groups.csv", groups_table(result))
        for (table, file) in ((:from_filing, "-groups-from-the-filing.csv"), (:from_open, "-groups-from-the-next-open.csv"))
            rows = stacked(result, table)
            isempty(rows) || save(name * file, rows)
        end
    end
    files
end

# FOMC releases

"""
An FOMC release: the title of the first of what came out at its time, which kinds (statement, projections, minutes,
other), when (UTC), and the minutes each of its bars spans (1, or 5 for a release captured after 30 days).
"""
struct Release
    accession::String
    title::String
    kinds::Vector{String}
    time::DateTime
    bar_minutes::Int
    complete::Bool
end

"""An instrument at a release: its Yahoo Finance symbol, the unit of its move (bp or %) and what the move is multiplied by."""
struct Contract
    symbol::String
    unit::String
    scale::Float64
end

struct FomcStudy
    folder::String
    title::String
    instruments::Vector{String}
    periods::Vector{Period}
    releases::Vector{Release}
    # By release and instrument
    contracts::Dict{Tuple{String, String}, Contract}
end

function read_fomc_study(folder, settings, periods)
    releases = Release[]
    path = joinpath(folder, "events.csv")
    if isfile(path)
        file = CSV.File(path; types = String)
        has_minutes = :bar_minutes in propertynames(file)
        for row in Tables.rows(file)
            ismissing(row.accepted) && continue
            minutes = has_minutes && !ismissing(row.bar_minutes) ? parse(Int, row.bar_minutes) : 1
            push!(releases, Release(row.accession, coalesce(row.name, ""), split(coalesce(row.form, "")), utc(row.accepted), minutes,
                coalesce(row.status, "") == "complete"))
        end
    end
    contracts = Dict{Tuple{String, String}, Contract}()
    path = joinpath(folder, "instruments.csv")
    if isfile(path)
        for row in Tables.rows(CSV.File(path; types = String))
            contracts[(row.accession, row.instrument)] = Contract(row.symbol, row.unit, parse(Float64, row.scale))
        end
    end
    FomcStudy(folder, get(settings, "title", basename(folder)), get(settings, "instruments", String[]), periods, releases, contracts)
end

"""
An instrument's move from the price `from` to `to`: in basis points of the implied rate (100 less the price) for the
federal funds and SOFR futures, times the contract's scale (Kuttner's for FF1); in basis points of the yield for a
yield; and in per cent (log) for a price.
"""
function instrument_move(instrument::AbstractString, contract::Contract, from, to)
    occursin(r"^(FF\d|SR1|SFR\d)$", instrument) && return -(to - from) * 100 * contract.scale
    contract.unit == "bp" && return (to - from) * 100
    100 * log(to / from)
end

"""
Each instrument's move at a release from the start of the period, a value a minute from `-before` to `after` (a bar
of five minutes holds its value for five); nothing for an instrument without the prices at the start.
"""
function analyse_release(study::FomcStudy, period::Period, release::Release, now::DateTime)
    prices = read_bars(joinpath(study.folder, "prices", release.accession * ".csv"), release.bar_minutes)
    known = release.complete ? typemax(DateTime) : now - SETTLE
    start = release.time - Minute(period.before)
    Dict(map(study.instruments) do instrument
        contract = get(study.contracts, (release.accession, instrument), nothing)
        bars = contract === nothing ? nothing : get(prices, contract.symbol, nothing)
        from = bars === nothing ? nothing : price_at(bars, start)
        from === nothing && return instrument => nothing
        function move_at(m)
            t = release.time + Minute(m)
            price = price_at(bars, t)
            t > known || price === nothing ? missing : instrument_move(instrument, contract, from, price)
        end
        instrument => Union{Missing, Float64}[move_at(m) for m in -period.before:period.after]
    end)
end

"""A period's results in a study of FOMC releases."""
struct FomcResult
    period::Period
    # A release a row: its accession, time, title and kinds, and a column per instrument with its move at the end of the period
    moves::NamedTuple
    # An instrument a row: its unit, how many releases have its move, and their mean, standard deviation and t
    summary::NamedTuple
    # By instrument, the mean move a minute at a time from the release, over the releases known at that minute, and its t
    paths::Dict{String, NamedTuple}
end

"""
    analyse(study::FomcStudy; now = Dates.now(UTC))

For each period of the study: each instrument's move at each release (`moves`, a table a release a row), the mean move
with its cross-sectional t (`summary`), and the mean move a minute at a time (`paths`, by instrument). Releases still
capturing prices count up to the last complete minute before `now`.
"""
function analyse(study::FomcStudy; now::DateTime = Dates.now(UTC))
    map(study.periods) do period
        minutes = -period.before:period.after
        paths = [analyse_release(study, period, release, now) for release in study.releases]
        last_move(path) = path === nothing ? missing : path[end]
        columns = [Symbol(instrument) => [last_move(p[instrument]) for p in paths] for instrument in study.instruments]
        moves = (; accession = [r.accession for r in study.releases], time = [r.time for r in study.releases],
            title = [r.title for r in study.releases], kinds = [join(r.kinds, " ") for r in study.releases], columns...)
        values = [known(last_move(p[instrument]) for p in paths) for instrument in study.instruments]
        unit(instrument) = something((c.unit for ((_, i), c) in study.contracts if i == instrument)..., "")
        summary = (instrument = copy(study.instruments), unit = unit.(study.instruments), releases = length.(values),
            mean = [isempty(v) ? missing : mean(v) for v in values], sd = [length(v) < 2 ? missing : std(v) for v in values],
            t = tstat.(values))
        averaged = Dict(instrument => average([p[instrument] for p in paths if p[instrument] !== nothing], minutes)
            for instrument in study.instruments)
        FomcResult(period, moves, summary, averaged)
    end
end

function report(io::IO, results::Vector{FomcResult})
    for result in results
        period, summary = result.period, result.summary
        println(io, "\n", period.name, ": ", period.before, " minutes before the release to ", period.after, " after; ",
            length(result.moves.accession), " releases")
        @printf(io, "%-8s %4s %9s %10s %10s %8s\n", "", "unit", "releases", "mean", "sd", "t")
        for k in eachindex(summary.instrument)
            @printf(io, "%-8s %4s %9d %10s %10s %8s\n", summary.instrument[k], summary.unit[k], summary.releases[k],
                format(summary.mean[k]), format(summary.sd[k]), format(summary.t[k]))
        end
    end
end

"""
    write_results(results::Vector{FomcResult}, folder)

Writes, for each period, the moves (`<period>-moves.csv`, a release a row), their means (`<period>-summary.csv`) and
the mean moves a minute at a time (`<period>-paths.csv`, an instrument and a minute a row) to `folder`; returns the files.
"""
function write_results(results::Vector{FomcResult}, folder::AbstractString)
    mkpath(folder)
    files = String[]
    save(name, table) = push!(files, CSV.write(joinpath(folder, name), table))
    for result in results
        name = slug(result.period.name)
        save(name * "-moves.csv", result.moves)
        save(name * "-summary.csv", result.summary)
        instruments = sort!(collect(keys(result.paths)))
        save(name * "-paths.csv", (instrument = [i for i in instruments for _ in result.paths[i].minute],
            minute = reduce(vcat, [result.paths[i].minute for i in instruments]; init = Int[]),
            releases = reduce(vcat, [result.paths[i].events for i in instruments]; init = Int[]),
            mean = reduce(vcat, [result.paths[i].mean for i in instruments]; init = Union{Missing, Float64}[]),
            t = reduce(vcat, [result.paths[i].t for i in instruments]; init = Union{Missing, Float64}[])))
    end
    files
end

include("bootstrap.jl")
include("panel.jl")
include("callaway.jl")

heading(study::Study) = string(study.title, ": benchmark ", study.benchmark, ", ",
    study.market_model ? "market model on the $(study.estimation_days) trading days before" : "market-adjusted", "; ",
    length(study.events), " events with a ticker and an acceptance time")
heading(study::FomcStudy) = string(study.title, ": FOMC releases; ", join(study.instruments, ", "), "; ", length(study.releases), " releases")

"""
    run_study(folder = pwd(); group_by = nothing, now = Dates.now(UTC))

Reads the study in `folder`, analyses it, prints the report and writes the results to `folder/results`; returns them.
A study of filings is grouped by `group_by` (see `analyse`), or else as its `study.toml` says.
"""
function run_study(folder::AbstractString = pwd(); group_by = nothing, now::DateTime = Dates.now(UTC))
    study = read_study(folder)
    println(heading(study))
    results = study isa Study ? analyse(study; group_by = something(group_by, study.group_by), now) : analyse(study; now)
    report(results)
    write_results(results, joinpath(folder, "results"))
    results
end

end
