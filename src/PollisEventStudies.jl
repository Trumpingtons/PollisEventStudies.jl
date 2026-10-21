"""
    PollisEventStudies

Intraday event studies of SEC filings, on the folders Pollis keeps for them: `study.toml` (the benchmark, the model of
normal returns and the periods), `events.csv` (the filings) and `prices/` (a file per filing with the minute bars of
the stock and of the benchmark). For each period: the average and cumulative average abnormal returns a minute at a
time with their cross-sectional t, the Boehmer, Musumeci and Poulsen (1991) test on standardized CARs, and abnormal
volume.

    results = run_study("~/.pollis/event-studies/my-study")
"""
module PollisEventStudies

using CSV, Tables, TOML, Dates, Statistics, Printf

export read_study, analyse, report, write_results, run_study

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
A filing: when EDGAR accepted it and the event's time (UTC), which is the earlier of the newswire release and the
acceptance when the study takes the first public time; whether the event came during the regular session, and the
next open after it.
"""
struct Event
    accession::String
    cik::String
    name::String
    ticker::String
    form::String
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
    periods = [Period(get(p, "name", "Period $k"), get(p, "before", 30), get(p, "after", 120), get(p, "next_open", false))
               for (k, p) in enumerate(get(settings, "periods", Dict{String, Any}[]))]
    isempty(periods) && push!(periods, Period("Default", 30, 120, true))
    events = Event[]
    path = joinpath(folder, "events.csv")
    if isfile(path)
        file = CSV.File(path; types = String)
        # An events.csv written before Pollis looked for newswire releases has no event_time: the acceptance is the event's time
        has_time = :event_time in propertynames(file)
        for row in Tables.rows(file)
            (ismissing(row.ticker) || ismissing(row.accepted)) && continue
            during = coalesce(row.during_session, "")
            time = has_time && !ismissing(row.event_time) ? row.event_time : row.accepted
            push!(events, Event(row.accession, coalesce(row.cik, ""), coalesce(row.name, ""), row.ticker, coalesce(row.form, ""),
                utc(row.accepted), utc(time), during == "true" ? true : during == "false" ? false : nothing,
                ismissing(row.next_open) ? nothing : utc(row.next_open), coalesce(row.status, "") == "complete"))
        end
    end
    Study(folder, get(settings, "title", basename(folder)), get(settings, "benchmark", "^GSPC"),
        get(settings, "model", "marketAdjusted") == "marketModel", get(settings, "estimation_days", 5), periods, events)
end

"""The minute bars of an event's prices file, by symbol."""
function read_prices(study::Study, event::Event)
    path = joinpath(study.folder, "prices", event.accession * ".csv")
    series = Dict{String, Series}()
    isfile(path) || return series
    rows = Tables.rowtable(CSV.File(path; types = Dict(:time => String, :session => String, :symbol => String)))
    for row in sort(rows; by = row -> row.time)
        bars = get!(() -> Series(DateTime[], Float64[], Float64[], Bool[]), series, row.symbol)
        push!(bars.ends, utc(row.time) + Minute(1))
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

"""How many of the bars' minutes end after from and at or before to."""
count_between(bars::Series, from::DateTime, to::DateTime) = max(0, searchsortedlast(bars.ends, to) - searchsortedlast(bars.ends, from))

"""The minute log returns of the stock and of the benchmark in the regular sessions of the `days` trading days before `start`."""
function estimation_returns(stock::Series, benchmark::Series, start::DateTime, days)
    # The benchmark trades only in the regular session: its minutes are the session's, all in one UTC day
    n = searchsortedlast(benchmark.ends, start - Millisecond(1))
    ends = benchmark.ends[1:n]
    chosen = Set(unique(Date.(ends))[max(1, end - days + 1):end])
    stock_returns, benchmark_returns = Float64[], Float64[]
    for k in 2:n
        (ends[k] - ends[k - 1] == Minute(1) && Date(ends[k]) in chosen) || continue
        price, before = price_at(stock, ends[k]), price_at(stock, ends[k - 1])
        (price === nothing || before === nothing) && continue
        push!(stock_returns, log(price / before))
        push!(benchmark_returns, log(benchmark.close[k] / benchmark.close[k - 1]))
    end
    stock_returns, benchmark_returns
end

"""Alpha (per minute of the regular session), beta and the residual standard deviation; nothing with too few returns."""
function fit(stock_returns, benchmark_returns)
    (length(stock_returns) < MIN_RETURNS || var(benchmark_returns) == 0) && return nothing
    beta = cov(benchmark_returns, stock_returns) / var(benchmark_returns)
    alpha = mean(stock_returns) - beta * mean(benchmark_returns)
    (alpha = alpha, beta = beta, sigma = std(stock_returns .- alpha .- beta .* benchmark_returns))
end

"""
An event's cumulative abnormal returns over a period, in per cent from the start of the period, a value a minute from
the filing and, for a filing outside the regular session in a period that follows the next open, from the next open;
and at the end of the period its CAR, standardized CAR, abnormal volume and beta. Nothing without the prices at its start.
"""
function analyse_event(study::Study, period::Period, event::Event, now::DateTime)
    prices = read_prices(study, event)
    stock, benchmark = get(prices, event.ticker, nothing), get(prices, study.benchmark, nothing)
    (stock === nothing || benchmark === nothing) && return nothing
    start = event.time - Minute(period.before)
    stock_start, benchmark_start = price_at(stock, start), price_at(benchmark, start)
    (stock_start === nothing || benchmark_start === nothing) && return nothing

    stock_returns, benchmark_returns = estimation_returns(stock, benchmark, start, study.estimation_days)
    model = study.market_model ? fit(stock_returns, benchmark_returns) : nothing
    sigma = model !== nothing ? model.sigma : length(stock_returns) >= MIN_RETURNS ? std(stock_returns .- benchmark_returns) : nothing
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

    finish = event.time + Minute(period.after)
    scar = sigma === nothing || ismissing(from_filing[end]) ? missing :
        from_filing[end] / 100 / (sigma * sqrt(max(1, count_between(stock, start, finish))))

    # The mean log(1 + volume) of the regular-session minutes after the event, or after the next open, less that of
    # the regular-session minutes before the event
    from = follows_open ? event.next_open : event.time
    window = [log1p(stock.volume[k]) for k in eachindex(stock.ends) if stock.regular[k] && from < stock.ends[k] <= min(from + Minute(period.after), known)]
    baseline = [log1p(stock.volume[k]) for k in eachindex(stock.ends) if stock.regular[k] && stock.ends[k] <= event.time]
    volume = isempty(window) || length(baseline) < MIN_VOLUME_MINUTES ? missing : mean(window) - mean(baseline)

    (accession = event.accession, ticker = event.ticker, name = event.name, from_filing = from_filing, from_open = from_open,
        car = from_filing[end], scar = scar, abnormal_volume = volume, beta = model === nothing ? missing : model.beta)
end

"""The cross-sectional t of a sample: its mean over its standard error; missing with fewer than two."""
tstat(values) = length(values) < 2 ? missing : mean(values) / (std(values) / sqrt(length(values)))

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
    analyse(study; now = Dates.now(UTC))

For each period of the study: its events' results (`events`, a table), each event's CAR a minute at a time from
the filing (`curves`), the AAR and CAAR from the filing
(`from_filing`, a table a minute a row) and from the next open (`from_open`, or nothing), and the tests at the end of
the period (`tests`). Events still capturing prices count up to the last complete minute before `now`.
"""
function analyse(study::Study; now::DateTime = Dates.now(UTC))
    map(study.periods) do period
        results = [result for result in (analyse_event(study, period, event, now) for event in study.events) if result !== nothing]
        opens = [result.from_open for result in results if result.from_open !== nothing]
        cars, scars, volumes = known(r.car for r in results), known(r.scar for r in results), known(r.abnormal_volume for r in results)
        mean_or_missing(values) = isempty(values) ? missing : mean(values)
        (period = period,
            events = (accession = [r.accession for r in results], ticker = [r.ticker for r in results], name = [r.name for r in results],
                car = [r.car for r in results], scar = [r.scar for r in results], abnormal_volume = [r.abnormal_volume for r in results],
                beta = [r.beta for r in results]),
            curves = [r.from_filing for r in results],
            from_filing = isempty(results) ? nothing : averages([r.from_filing for r in results], -period.before:period.after),
            from_open = isempty(opens) ? nothing : averages(opens, 0:period.after),
            open_events = length(opens),
            tests = (car = mean_or_missing(cars), car_t = tstat(cars), car_events = length(cars),
                bmp = mean_or_missing(scars), bmp_t = tstat(scars), bmp_events = length(scars),
                abnormal_volume = mean_or_missing(volumes), abnormal_volume_t = tstat(volumes), abnormal_volume_events = length(volumes)))
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
        @printf(io, "\nAt +%d minutes: CAR %s%% (t = %s, %d events)\n", period.after, format(tests.car), format(tests.car_t), tests.car_events)
        if tests.bmp_events == 0
            println(io, "Boehmer, Musumeci and Poulsen: needs the estimation days, which only the market model captures")
        else
            @printf(io, "Boehmer, Musumeci and Poulsen: mean standardized CAR %s (t = %s, %d events)\n", format(tests.bmp), format(tests.bmp_t), tests.bmp_events)
        end
        @printf(io, "Abnormal volume: %s (t = %s, %d events)\n", format(tests.abnormal_volume), format(tests.abnormal_volume_t), tests.abnormal_volume_events)
    end
end

slug(name) = strip(replace(lowercase(name), r"[^a-z0-9]+" => "-"), '-')

"""
    write_results(results, folder)

Writes, for each period, its events (`<period>-events.csv`) and its AAR and CAAR a minute a row
(`<period>-from-the-filing.csv`, `<period>-from-the-next-open.csv`) to `folder`; returns the files.
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
    end
    files
end

"""
    run_study(folder = pwd(); now = Dates.now(UTC))

Reads the study in `folder`, analyses it, prints the report and writes the results to `folder/results`; returns them.
"""
function run_study(folder::AbstractString = pwd(); now::DateTime = Dates.now(UTC))
    study = read_study(folder)
    println(study.title, ": benchmark ", study.benchmark, ", ",
        study.market_model ? "market model on the $(study.estimation_days) trading days before" : "market-adjusted", "; ",
        length(study.events), " events with a ticker and an acceptance time")
    results = analyse(study; now)
    report(results)
    write_results(results, joinpath(folder, "results"))
    results
end

end
