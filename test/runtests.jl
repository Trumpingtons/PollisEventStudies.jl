using PollisEventStudies, Test, Random, Dates, CSV, Tables, WildBootTests, MixedModels

const BENCHMARK = "^GSPC"
const DAYS = Date.(2026, 10, [12, 13, 14, 15, 16, 19, 20])
const SESSION_MINUTES = 390

iso(t::DateTime) = Dates.format(t, dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
open_of(day::Date) = DateTime(day) + Hour(14) + Minute(30)

"""
A study folder with synthetic minute prices: the stock's return is alpha plus beta times the benchmark's plus noise,
with a jump at the minute EDGAR accepted the filing (or at the next open), and five times the volume in the hour after.
"""
function synthetic_study(folder; model = "marketModel", events, status = "complete")
    mkpath(joinpath(folder, "prices"))
    write(joinpath(folder, "study.toml"), """
        title = "Synthetic"
        alert = "test"
        benchmark = "$BENCHMARK"
        model = "$model"
        estimation_days = 5

        [[periods]]
        name = "Default"
        before = 30
        after = 120
        next_open = true

        [[periods]]
        name = "Short Run"
        before = 5
        after = 15
        next_open = false
        """)
    rng = Xoshiro(1)
    open(joinpath(folder, "events.csv"), "w") do io
        println(io, "accession,cik,name,ticker,form,items,filed,accepted,during_session,next_open,status,url")
        for e in events
            next_open = e.next_open === nothing ? "" : iso(e.next_open)
            println(io, "$(e.accession),1,\"Company $(e.ticker), Inc.\",$(e.ticker),$(get(e, :form, "8-K")),$(get(e, :items, "2.02")),2026-10-19,$(iso(e.accepted)),$(e.during),$next_open,$status,https://www.sec.gov/")
        end
    end
    for e in events
        open(joinpath(folder, "prices", e.accession * ".csv"), "w") do io
            println(io, "time,session,symbol,open,high,low,close,volume")
            stock, level = 100.0, 5000.0
            jumps(t) = (t <= e.accepted < t + Minute(1)) || t == e.next_open
            hot(t) = Minute(0) <= t - e.accepted < Minute(60) || (e.next_open !== nothing && Minute(0) <= t - e.next_open < Minute(60))
            for day in DAYS, m in 0:SESSION_MINUTES - 1
                t = open_of(day) + Minute(m)
                market = 0.0005 * randn(rng)
                level *= exp(market)
                stock *= exp(0.00001 + e.beta * market + 0.001 * randn(rng) + (jumps(t) ? e.jump : 0.0))
                volume = rand(rng, 1000:2000) * (hot(t) ? 5 : 1)
                # A thinly traded stock has a bar only in the minutes it trades in
                m % get(e, :trades_every, 1) == 0 && println(io, "$(iso(t)),regular,$(e.ticker),$stock,$stock,$stock,$stock,$volume")
                println(io, "$(iso(t)),regular,$BENCHMARK,$level,$level,$level,$level,0")
            end
        end
    end
    folder
end

"""
A study folder of FOMC releases with synthetic bars: each instrument's price is flat but for a jump in the bar that
starts at the release; `releases` give the jumps, the bars' minutes and the status.
"""
function synthetic_fomc_study(folder; releases)
    mkpath(joinpath(folder, "prices"))
    write(joinpath(folder, "study.toml"), """
        title = "FOMC"
        source = "fomc"
        releases = ["statement", "projections"]
        instruments = ["FF1", "SFR1", "FVX", "^GSPC"]

        [[periods]]
        name = "Tight"
        before = 10
        after = 20
        """)
    open(joinpath(folder, "events.csv"), "w") do io
        println(io, "accession,cik,name,ticker,form,items,filed,accepted,during_session,next_open,status,url,event_time,bar_minutes")
        for r in releases
            println(io, "$(r.accession),,Federal Reserve issues FOMC statement,,statement projections,,$(Date(r.time)),$(iso(r.time)),,,$(r.status),https://www.federalreserve.gov/,$(iso(r.time)),$(r.minutes == 5 ? 5 : "")")
        end
    end
    open(joinpath(folder, "instruments.csv"), "w") do io
        println(io, "accession,event_time,instrument,symbol,unit,scale,bar_minutes")
        for r in releases
            println(io, "$(r.accession),$(iso(r.time)),FF1,ZQ$(r.accession).CBT,bp,$(r.scale),$(r.minutes)")
            println(io, "$(r.accession),$(iso(r.time)),SFR1,SR3$(r.accession).CME,bp,1,$(r.minutes)")
            println(io, "$(r.accession),$(iso(r.time)),FVX,^FVX,bp,1,$(r.minutes)")
            println(io, "$(r.accession),$(iso(r.time)),^GSPC,^GSPC,%,1,$(r.minutes)")
        end
    end
    for r in releases
        open(joinpath(folder, "prices", r.accession * ".csv"), "w") do io
            println(io, "time,session,symbol,open,high,low,close,volume")
            # No bars for SFR1: its contract is gone from Yahoo Finance
            for (symbol, before, after) in (("ZQ$(r.accession).CBT", r.ff...), ("^FVX", r.fvx...), ("^GSPC", r.gspc...))
                for m in -60:r.minutes:55
                    t = r.time + Minute(m)
                    price = m < 0 ? before : after
                    println(io, "$(iso(t)),regular,$symbol,$price,$price,$price,$price,0")
                end
            end
        end
    end
    folder
end

const RELEASES = [
    # FF1 down 10 cents, scaled by 30 / (30 - 10) for a release on the 10th: 15 bp; the 5-year yield up 5 bp; the S&P 500 down 1%
    (accession = "fomc-20260910T1800", time = DateTime(2026, 9, 10, 18), minutes = 1, status = "complete", scale = 1.5,
        ff = (95.0, 94.9), fvx = (4.0, 4.05), gspc = (5000.0, 5000 * exp(-0.01))),
    # Older, in five-minute bars: FF1 down 5 cents, 5 bp; the yield down 3 bp; the S&P 500 up 0.5%
    (accession = "fomc-20260729T1800", time = DateTime(2026, 7, 29, 18), minutes = 5, status = "complete", scale = 1.0,
        ff = (95.5, 95.45), fvx = (4.2, 4.17), gspc = (5000.0, 5000 * exp(0.005))),
]

const EVENTS = [
    (accession = "E1", ticker = "AAA", accepted = DateTime(2026, 10, 19, 16, 0, 30), during = "true", next_open = nothing, jump = 0.02, beta = 1.2),
    # After the close on Friday: the jump comes at Monday's open
    (accession = "E2", ticker = "BBB", accepted = DateTime(2026, 10, 16, 21, 30, 10), during = "false", next_open = DateTime(2026, 10, 19, 14, 30), jump = -0.03, beta = 0.8),
    (accession = "E3", ticker = "CCC", accepted = DateTime(2026, 10, 20, 15, 0, 0), during = "true", next_open = nothing, jump = 0.01, beta = 1.0),
]

"""
A panel of `units` units over the years 1 to 12: `cohort(i)` the year unit i is first treated (0 for never) and
`effect(e, l)` the effect on cohort e at event time l, from 0, or before it with a pre-trend; the noise is a function
of the unit and the year, the same on every machine.
"""
function synthetic_panel(cohort, effect; units = 60, years = 1:12)
    rows = (; unit = String[], year = Int[], y = Float64[], first = Int[], w = Float64[], x = Float64[], region = String[])
    for i in 1:units, t in years
        e = cohort(i)
        x = cos(1.3i + 0.7t)
        push!(rows.unit, "u$i"); push!(rows.year, t); push!(rows.first, e); push!(rows.x, x)
        push!(rows.w, 1 + 0.5 * mod(i, 3)); push!(rows.region, "r$(mod(i, 12))")
        push!(rows.y, sin(i) + 0.2t + 0.5x + (e == 0 ? 0.0 : effect(e, t - e)) + 0.05 * sin(3.1i * t + i))
    end
    rows
end
# Three cohorts and units never treated, with effects that differ from cohort to cohort
STAGGERED = synthetic_panel(i -> i <= 15 ? 4 : i <= 30 ? 7 : i <= 45 ? 9 : 0,
    (e, l) -> l < 0 ? 0.0 : e == 4 ? 0.4 * (l + 1) : e == 7 ? 2.0 : -1.0)
# The same, unbalanced: the even units miss a year or two
GAPS = map(c -> c[[(i = parse(Int, STAGGERED.unit[k][2:end]); isodd(i) || mod(7i + 3 * STAGGERED.year[k], 11) != 0)
    for k in eachindex(STAGGERED.unit)]], STAGGERED)

@testset "PollisEventStudies" begin
    @testset "market model" begin
        folder = synthetic_study(mktempdir(); events = EVENTS)
        study = read_study(folder)
        @test (study.title, study.benchmark, study.market_model, length(study.periods), length(study.events)) == ("Synthetic", BENCHMARK, true, 2, 3)

        results = analyse(study)
        default = results[1]
        # The betas are recovered, and each jump shows in its event's CAR the minute after
        @test all(abs.(default.events.beta .- [1.2, 0.8, 1.0]) .< 0.15)
        minute_one = 30 + 2
        @test abs(default.curves[1][minute_one] - 2.0) < 0.5
        @test abs(default.curves[3][minute_one] - 1.0) < 0.5
        # The filing after the close is followed from the next open, where its jump is
        @test default.open_events == 1
        @test default.from_open.caar[2] - default.from_open.caar[1] < -2.5
        @test default.from_filing.events[1] == 3 && default.from_filing.caar[1] == 0
        # Its CAR at the end of the period is that of the period's minutes after the next open, not of the closed market
        # after the filing; in a period that does not follow the next open, the latter
        @test default.events.car[2] == default.from_open.caar[end] < -1 && results[2].events.car[2] == 0
        @test default.tests.bmp_events == 3 && default.tests.abnormal_volume_events == 3
        # The filing after Friday's close is followed to Monday afternoon, over Monday's filing: their abnormal returns,
        # independent, are correlated over their estimation days in common by a little, and the adjusted t moves by as little
        @test 0 < abs(default.tests.correlation) < 0.05 && abs(default.tests.bmp_kp_t - default.tests.bmp_t) < 0.05
        @test isfinite(default.tests.grank_t) && 0 <= default.tests.grank_p <= 1
        # Five times the volume for the first 60 of the 120 minutes: log(5) / 2
        @test abs(default.tests.abnormal_volume - log(5) / 2) < 0.1
        # The jump in the minute after the filing raises the abnormal returns' realized volatility over the next 120 minutes
        @test default.events.abnormal_volatility[1] > 0.4 && default.tests.abnormal_volatility_events == 3
        @test default.tests.gsign_events == 3 && 0 <= default.tests.gsign_p <= 1 && 0 <= default.tests.wilcoxon_p <= 1
        # A period that does not follow the next open has no curve from it
        @test results[2].from_open === nothing

        files = write_results(results, joinpath(folder, "results"))
        @test basename.(files) == ["default-events.csv", "default-outcomes.csv", "default-from-the-filing.csv", "default-from-the-next-open.csv",
            "short-run-events.csv", "short-run-outcomes.csv", "short-run-from-the-filing.csv"]
        @test length(CSV.File(files[3])) == 30 + 120 + 1

        output = sprint(report, results)
        @test occursin("Minutes from the next open, 1 events", output)
        @test occursin("Boehmer, Musumeci and Poulsen: mean standardized CAR", output)
    end

    @testset "market-adjusted" begin
        study = read_study(synthetic_study(mktempdir(); model = "marketAdjusted", events = EVENTS))
        default = analyse(study)[1]
        @test all(ismissing, default.events.beta)
        # The prices before the period are those of the market model's estimation days: here they are there anyway
        @test default.tests.bmp_events == 3
        @test abs(default.curves[1][32] - 2.0) < 0.5
    end

    @testset "thin trading" begin
        # A stock that trades every fourth minute: with a return of zero in the minutes between, its beta comes out at
        # 0.33 here; measured from trade to trade, 0.93
        thin = merge(EVENTS[3], (trades_every = 4,))
        default = analyse(read_study(synthetic_study(mktempdir(); events = [thin])))[1]
        @test abs(default.events.beta[1] - 1.0) < 0.25
        @test default.tests.bmp_events == 1
    end

    @testset "sign and signed-rank tests" begin
        P = PollisEventStudies
        # Two of four estimation returns positive: half the CARs are expected positive, and three of four are
        cowan = P.generalized_sign([(car = car, sars = [1.0, -1.0, 2.0, -2.0]) for car in [1.0, 2.0, 3.0, -1.0]])
        @test cowan == (z = 1.0, positive = 3)
        # With one event, no z, as the cross-sectional t and EventStudy's GeneralizedSignTest
        one = P.generalized_sign([(car = 1.0, sars = [1.0, -1.0])])
        @test ismissing(one.z) && one.positive == 1
        # R's wilcox.test(x, y, paired = TRUE) of its own example, exact: V = 40, p = 0.03906
        x = [1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30]
        y = [0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29]
        exact = P.signed_rank(x .- y)
        @test exact.v == 40 && round(exact.p; digits = 5) == 0.03906
        # Sixty values, two of them zero and the others tied in pairs: R's normal approximation, V = 1340, p = 0.0001777228713
        z = [-1; repeat(-0.9:0.1:1.9; inner = 2); 2]
        tied = P.signed_rank(z)
        @test tied.v == 1340 && round(tied.p; sigdigits = 8) == round(0.0001777228713; sigdigits = 8)
    end

    @testset "outcomes" begin
        day = DateTime(DAYS[1])
        events = [
            # Its prices file has no bars for the ticker in events.csv
            merge(EVENTS[1], (accession = "E4",)),
            # Ten minutes after the first open: the period starts before the first bar
            merge(EVENTS[3], (accession = "E5", accepted = day + Hour(14) + Minute(40))),
            # 45 minutes after the first open: the 14 minutes before the period are too few for the market model
            merge(EVENTS[3], (accession = "E6", accepted = day + Hour(15) + Minute(15))),
            EVENTS[3]]
        folder = synthetic_study(mktempdir(); events)
        path = joinpath(folder, "events.csv")
        write(path, replace(read(path, String), "E4,1,\"Company AAA, Inc.\",AAA" => "E4,1,\"Company AAA, Inc.\",ZZZ"))
        default = analyse(read_study(folder))[1]
        @test (default.outcomes.accession, default.outcomes.outcome, default.outcomes.notes) == (["E4", "E5", "E6", "E3"],
            ["no_stock_prices", "no_start_price", "used", "used"], ["", "", "market_adjusted no_scar", ""])
        @test default.events.accession == ["E6", "E3"] && default.tests.bmp_events == 1
        output = sprint(report, [default])
        @test occursin("2 events with prices\n  left out, no price at the start of the period: 1\n  left out, no prices for the stock: 1\n" *
            "  market-adjusted: too few estimation minutes to fit the market model: 1\n", output)
    end

    @testset "Kolari and Pynnonen" begin
        P = PollisEventStudies
        rng = Xoshiro(2)
        times = DateTime(2026, 10, 19, 14, 31) .+ Minute.(0:199)
        shared = randn(rng, 200)
        window(hour) = (DateTime(2026, 10, 20, hour), DateTime(2026, 10, 20, hour + 2))
        event(sars, scar, hour) = (sar_times = times, sars = sars, scar = scar, window = window(hour))
        # Two events at the same time with the same abnormal returns, and one later whose period does not overlap
        results = [event(shared, 1.0, 14), event(shared, 1.0, 15), event(randn(rng, 200), 1.0, 18)]
        @test P.cross_correlation(results) == 1 / 3
        @test P.correlation_factor(1 / 3, 3) == sqrt((2 / 3) / (5 / 3))
        @test P.correlation_factor(0.0, 3) == 1

        # Each event's CAR far above its estimation days' abnormal returns: the highest rank, and a t well above 3
        far = [event(randn(rng, 200), 10.0 + k, 14 + 2k) for k in 0:4]
        grank = P.generalized_rank(far)
        @test grank.df == 201 - 2 && grank.t > 3 && P.pvalue(grank.t, grank.df) < 0.01
        # Among them: no evidence
        near = [event(randn(rng, 200), 0.1 * (-1)^k, 14 + 2k) for k in 0:4]
        @test abs(P.generalized_rank(near).t) < 2
        @test ismissing(P.generalized_rank(far[1:1]).t)
    end

    @testset "groups" begin
        events = [merge(EVENTS[1], (form = "8-K", items = "2.02")), merge(EVENTS[2], (form = "10-Q", items = "")),
            merge(EVENTS[3], (form = "8-K", items = "2.02 9.01"))]
        folder = synthetic_study(mktempdir(); events)
        study = read_study(folder)
        @test study.group_by == "" && study.events[3].items == ["2.02", "9.01"]
        @test isempty(analyse(study)[1].groups)

        default = analyse(study; group_by = "item")[1]
        @test [(g.name, g.events.accession) for g in default.groups] == [("Item 2.02", ["E1", "E3"]), ("Item 9.01", ["E3"]), ("No item", ["E2"])]
        # A group's tests are those of its events
        @test default.groups[3].tests.car == default.events.car[2]
        @test [g.name for g in analyse(study; group_by = "session")[1].groups] == ["During the session", "Outside the session"]
        @test [g.name for g in analyse(study; group_by = event -> event.ticker == "AAA" ? "A" : "Others")[1].groups] == ["A", "Others"]
        @test_throws ArgumentError analyse(study; group_by = "colour")

        results = analyse(study; group_by = "form")
        @test [length(g.events.accession) for g in results[1].groups] == [1, 2]
        output = sprint(report, results)
        @test occursin("By form, at +120 minutes", output) && occursin("10-Q", output)
        files = basename.(write_results(results, joinpath(folder, "results")))
        @test "default-groups.csv" in files && "default-groups-from-the-filing.csv" in files && "default-groups-from-the-next-open.csv" in files
        @test length(CSV.File(joinpath(folder, "results", "default-groups-from-the-filing.csv"))) == 2 * (30 + 120 + 1)
    end

    @testset "events still capturing" begin
        study = read_study(synthetic_study(mktempdir(); events = EVENTS[1:1], status = "capturing"))
        # Ten minutes after the filing, the CARs are known to the last complete minute and no further
        default = analyse(study; now = EVENTS[1].accepted + Minute(10))[1]
        known = findall(!ismissing, default.curves[1])
        @test last(known) == 30 + 1 + 8
        @test default.from_filing.events[end] == 0
        @test ismissing(default.tests.car)
    end

    @testset "FOMC releases" begin
        folder = synthetic_fomc_study(mktempdir(); releases = RELEASES)
        study = read_study(folder)
        @test (study.title, study.instruments, length(study.periods), length(study.releases)) == ("FOMC", ["FF1", "SFR1", "FVX", "^GSPC"], 1, 2)
        @test study.releases[2].bar_minutes == 5 && study.releases[1].kinds == ["statement", "projections"]

        tight = analyse(study)[1]
        moves = tight.moves
        @test isapprox(collect(skipmissing(moves.FF1)), [15, 5]; atol = 1e-9)
        @test isapprox(collect(skipmissing(moves.FVX)), [5, -3]; atol = 1e-9)
        @test isapprox(collect(skipmissing(moves[Symbol("^GSPC")])), [-1, 0.5]; atol = 1e-9)
        # An instrument whose contract has no prices has no move
        @test all(ismissing, moves.SFR1)
        @test tight.summary.releases == [2, 0, 2, 2] && tight.summary.unit == ["bp", "bp", "bp", "%"]
        @test isapprox(tight.summary.mean[1], 10; atol = 1e-9) && ismissing(tight.summary.mean[2])

        # A minute after the release, the one-minute bars show the jump and the five-minute bars do not until five minutes after
        ff = tight.paths["FF1"]
        at(m) = findfirst(==(m), ff.minute)
        @test ff.events[at(0)] == 2 && isapprox(ff.mean[at(0)], 0; atol = 1e-9)
        @test isapprox(ff.mean[at(1)], 7.5; atol = 1e-9)
        @test isapprox(ff.mean[at(5)], 10; atol = 1e-9)
        @test tight.paths["SFR1"].events[at(0)] == 0

        files = write_results(analyse(study), joinpath(folder, "results"))
        @test basename.(files) == ["tight-moves.csv", "tight-summary.csv", "tight-paths.csv"]
        @test length(CSV.File(files[3])) == 4 * (10 + 20 + 1)
        output = sprint(report, analyse(study))
        @test occursin("Tight: 10 minutes before the release to 20 after; 2 releases", output)
    end

    @testset "FOMC release still capturing" begin
        release = merge(RELEASES[1], (status = "capturing",))
        study = read_study(synthetic_fomc_study(mktempdir(); releases = [release]))
        # Five minutes after the release, the moves are known to the last complete minute and no further
        tight = analyse(study; now = release.time + Minute(5))[1]
        known = findall(!ismissing, tight.paths["FVX"].mean)
        @test tight.paths["FVX"].minute[last(known)] == 3
        @test ismissing(tight.moves.FVX[1])
    end

    @testset "no events" begin
        folder = mktempdir()
        write(joinpath(folder, "study.toml"), "title = \"Empty\"\nalert = \"x\"\n")
        write(joinpath(folder, "events.csv"), "accession,cik,name,ticker,form,items,filed,accepted,during_session,next_open,status,url\n")
        results = analyse(read_study(folder))
        @test length(results) == 1 && results[1].period.name == "Default" && results[1].from_filing === nothing
        @test occursin("0 events with prices", sprint(report, results))
    end

    @testset "panel event studies" begin
        options = (; unit = "unit", time = "year", outcome = "y", cohort = "first", leads = 8, lags = 8)
        study(table; changes...) = panel_event_study(table; merge(options, changes)...)
        rounded(x) = round(x; digits = 6)
        cell(table, keys...) = findfirst(k -> all(getproperty(table, name)[k] == value for (name, value) in keys), eachindex(table.event_time))

        # R's fixest, feols(y ~ x + sunab(cohort, year, ref.p = -1) | unit + year, weights = ~w, cluster = ~region) and
        # its i(rel) two-way fixed effects, on the same panel: each cohort's effect, the TWFE effects and their standard errors
        weighted = study(STAGGERED; weights = "w", covariates = ["x"], cluster = "region")
        cohort_rows = [cell(weighted.cohorts, :cohort => c, :event_time => l) for (c, l) in ((9, -8), (4, -3), (7, -2), (4, 0), (7, 0))]
        event_rows = [cell(weighted.effects, :event_time => l) for l in (-8, -2, 0, 4, 8)]
        @test [rounded.((weighted.cohorts.effect[k], weighted.cohorts.se[k])) for k in cohort_rows] ==
            [(0.010423, 0.019092), (0.008481, 0.011561), (-0.013831, 0.01562), (0.395517, 0.019796), (1.990594, 0.013136)]
        @test [rounded.((weighted.effects.twfe_effect[k], weighted.effects.twfe_se[k])) for k in event_rows] ==
            [(1.316918, 0.036649), (0.03789, 0.015169), (0.456324, 0.056954), (1.759524, 0.047242), (4.162545, 0.050693)]
        # and on the panel with gaps
        gappy = study(GAPS; weights = "w", covariates = ["x"], cluster = "region")
        @test ([(k = cell(gappy.effects, :event_time => l); rounded.((gappy.effects.twfe_effect[k], gappy.effects.twfe_se[k])))
                for l in (-8, 0, 8)],
            [(k = cell(gappy.cohorts, :cohort => c, :event_time => l); rounded.((gappy.cohorts.effect[k], gappy.cohorts.se[k])))
                for (c, l) in ((9, -8), (7, 0))]) ==
            ([(1.305337, 0.037126), (0.449751, 0.055236), (4.139462, 0.047372)], [(0.004666, 0.021067), (1.985342, 0.012569)])
        # and its aggregated effects, the cohorts' shares held fixed
        fixed = study(STAGGERED; weights = "w", covariates = ["x"], cluster = "region", fixed_shares = true)
        rows = [cell(fixed.effects, :event_time => l) for l in (-2, 0, 3, 8)]
        @test [rounded.((fixed.effects.sa_effect[k], fixed.effects.sa_se[k])) for k in rows] ==
            [(-0.00045, 0.012806), (0.464758, 0.010749), (0.862334, 0.015045), (3.608945, 0.020174)]

        # R's did, att_gt(control_group = "nevertreated", base_period = "universal") and aggte(type = "dynamic"),
        # without a small-sample factor: the effect at each event time, and its standard error with the variance of the
        # cohorts' shares
        unweighted = study(STAGGERED; small_sample = false)
        rows = [cell(unweighted.effects, :event_time => l) for l in (-8, -2, 0, 3, 8)]
        @test [rounded.((unweighted.effects.sa_effect[k], unweighted.effects.sa_se[k])) for k in rows] ==
            [(-0.002116, 0.163785), (0.004404, 0.040237), (0.468011, 0.186927), (0.871297, 0.228423), (3.617982, 0.018838)]
        # and, with the normal as did, its p-values and 95% bands (aggte(bstrap = FALSE, cband = FALSE))
        normal = study(STAGGERED; small_sample = false, distribution = "normal")
        rows = [cell(normal.effects, :event_time => l) for l in (-2, 0, 3)]
        @test [rounded.((normal.effects.sa_p[k], normal.effects.sa_low[k], normal.effects.sa_high[k])) for k in rows] ==
            [(0.912849, -0.074459, 0.083267), (0.01229, 0.101641, 0.834380), (0.000137, 0.423597, 1.318997)]
        @test normal.pretrends.test == ["chi2", "chi2"] && all(ismissing, normal.pretrends.df2) && normal.pretrends.p[2] > 0.5
        @test_throws ArgumentError study(STAGGERED; distribution = "z")

        # Sup-t bands, simultaneous over the 16 event times: wider than pointwise, narrower than Bonferroni's (2.955);
        # 2.775 and 2.777 with 200,000 draws (did's multiplier bootstrap, its weights +1 or -1 over the 60 clusters, gives 2.708)
        uniform = study(STAGGERED; small_sample = false, distribution = "normal", uniform = true, draws = 20000, seed = 1)
        c = uniform.sample.critical
        @test isapprox(c.twfe, 2.775; atol = 0.03) && isapprox(c.sa, 2.777; atol = 0.03)
        known = findall(!ismissing, uniform.effects.sa_effect)
        @test isapprox(uniform.effects.sa_high[known] .- uniform.effects.sa_effect[known], c.sa .* uniform.effects.sa_se[known])
        @test isequal(uniform.effects.sa_effect, normal.effects.sa_effect)
        @test study(STAGGERED; small_sample = false, distribution = "normal", uniform = true, draws = 20000, seed = 1).sample.critical == c
        # Under t, heavier tails; with a single event time, the pointwise value
        @test study(STAGGERED; uniform = true, seed = 1).sample.critical.sa > c.sa
        @test study(STAGGERED; leads = 1, lags = 0, distribution = "normal", uniform = true).sample.critical.sa |> x -> round(x; digits = 6) == 1.959964
        @test_throws ArgumentError study(STAGGERED; uniform = true, draws = 10)
        # did's multiplier bootstrap in place of the sup-t: 2.706 with 200,000 draws, did's 2.708
        did_bands = study(STAGGERED; small_sample = false, distribution = "normal", uniform = true, multiplier = true,
            draws = 20000, seed = 1)
        @test isapprox(did_bands.sample.critical.sa, 2.708; atol = 0.03)
        @test_throws ArgumentError study(STAGGERED; multiplier = true)

        # R's fwildclusterboot, boottest(type = "rademacher", B = 9999) on fixest's i(rel) and sunab(no_agg = TRUE) with
        # the unit and year effects as indicators, Sun and Abraham's effect the shares' combination of the cohorts': with 12
        # clusters, all 4,096 patterns of signs, so the p-values are R's. The intervals' bounds are where the p-value, a step
        # function, falls below 0.05: TWFE's at 3, 1.158875; WildBootTests stops within a step of 0.05 (1.159517), R's
        # uniroot after 10 iterations (1.158985), both on the step past the bound
        wild = study(STAGGERED; weights = "w", covariates = ["x"], cluster = "region", fixed_shares = true, wild = true)
        rows = [cell(wild.effects, :event_time => l) for l in (-2, 0, 3)]
        @test wild.sample.wild_draws == 4096
        @test [rounded.((wild.effects.twfe_wild_p[k], wild.effects.sa_wild_p[k])) for k in rows] ==
            [(0.032715, 0.966309), (0.0, 0.0), (0.0, 0.0)]
        intervals = [(wild.effects.twfe_wild_low[k], wild.effects.twfe_wild_high[k], wild.effects.sa_wild_low[k],
            wild.effects.sa_wild_high[k]) for k in rows]
        @test all(isapprox.(collect.(intervals), [[0.003463, 0.072763, -0.030166, 0.026207],
            [0.326269, 0.586258, 0.441556, 0.488692], [0.916298, 1.158985, 0.828199, 0.892420]]; atol = 1e-3))
        @test wild.pretrends.wild_p[1] < 0.001 && wild.pretrends.wild_p[2] > 0.5
        @test_throws ArgumentError study(STAGGERED; wild = true, wild_weights = "mammen")
        # Webb's weights, against boottest(type = "webb", B = 999999): 6^12 patterns are too many to enumerate, so both
        # draw at random, and TWFE's p-values and intervals agree within the noise of the draws
        webb = study(STAGGERED; weights = "w", covariates = ["x"], cluster = "region", fixed_shares = true, wild = true,
            wild_weights = "webb", draws = 99999, seed = 1)
        rows = [cell(webb.effects, :event_time => l) for l in (-2, 0, 3)]
        @test webb.sample.wild_draws == 99999 && all(isapprox.(
            [[webb.effects.twfe_wild_p[k], webb.effects.twfe_wild_low[k], webb.effects.twfe_wild_high[k]] for k in rows],
            [[0.030189, 0.003807, 0.072266], [0.000086, 0.326607, 0.58471], [0.000017, 0.918462, 1.159012]]; atol = 2e-3))

        # Effects that differ from cohort to cohort: Sun and Abraham recovers them, averaged over the cohorts at each event
        # time; TWFE mixes them in, before treatment too
        result = study(STAGGERED)
        effect(c, l) = c == 4 ? 0.4 * (l + 1) : c == 7 ? 2.0 : -1.0
        present(l) = [c for c in (4, 7, 9) if 1 <= c + l <= 12]
        truth = Dict(l => l < 0 ? 0.0 : sum(effect(c, l) for c in present(l)) / length(present(l)) for l in result.effects.event_time)
        known = findall(!ismissing, result.effects.sa_effect)
        @test maximum(abs(result.effects.sa_effect[k] - truth[result.effects.event_time[k]]) for k in known) < 0.05
        @test result.pretrends.p[1] < 1e-6 && result.pretrends.p[2] > 0.5
        @test (result.sample.observations, result.sample.units, result.sample.cohorts, result.sample.never_treated) == (720, 60, [4, 7, 9], 15)

        # The same with a treatment column in place of the cohorts
        treated = merge(STAGGERED, (; d = [e != 0 && t >= e ? 1 : 0 for (e, t) in zip(STAGGERED.first, STAGGERED.year)]))
        @test isequal(study(treated; cohort = nothing, treatment = "d").effects, result.effects)

        # Without units never treated, those treated last are the control, before their treatment
        some = findall(!=(0), STAGGERED.first)
        last_treated = study(map(column -> column[some], STAGGERED))
        @test last_treated.sample.control == "treated last (9), before their treatment" && last_treated.sample.observations == 45 * 8
        @test_throws ArgumentError study(map(column -> column[some], STAGGERED); control = "never")

        # Units treated from the first year are left out; event times beyond the ends count at the ends
        early = merge(STAGGERED, (; first = [u in ("u1", "u2") ? 1 : e for (u, e) in zip(STAGGERED.unit, STAGGERED.first)]))
        binned = study(early; leads = 3, lags = 2)
        @test binned.sample.dropped == 2 && binned.effects.event_time == [-3, -2, 0, 1, 2]
        # At -3, the years 3 or more before treatment: a year of each of cohort 4's 13 units, 4 of cohort 7's 15
        # and 6 of cohort 9's; at 2, the years 2 or more after: 7, 4 and 2
        @test binned.effects.observations == [13 + 4 * 15 + 6 * 15, 43, 43, 43, 7 * 13 + 4 * 15 + 2 * 15]

        # A cohort whose outcome moves before its treatment fails Sun and Abraham's pre-trends test
        pretrend = synthetic_panel(i -> i <= 30 ? 6 : 0, (e, l) -> l < 0 ? 0.3 * l : 1.0)
        @test study(pretrend).pretrends.p[2] < 1e-6

        # A study folder
        folder = mktempdir()
        CSV.write(joinpath(folder, "panel.csv"), STAGGERED)
        write(joinpath(folder, "study.toml"), """
            title = "Staggered"
            source = "panel"
            unit = "unit"
            time = "year"
            outcome = "y"
            cohort = "first"
            leads = 8
            lags = 8
            estimators = ["panel_event_study", "callaway_santanna"]

            [callaway_santanna]
            control = "notyet"
            """)
        output = mktemp() do path, io
            redirect_stdout(() -> run_study(folder), io)
            close(io)
            read(path, String)
        end
        @test occursin("720 observations of 60 units in 60 clusters; cohorts 4, 7, 9; 15 never treated", output)
        @test occursin("Callaway and Sant'Anna: 60 units in 60 clusters, times 1 to 12; cohorts 4, 7, 9; control: not yet treated", output)
        @test sort(readdir(joinpath(folder, "results"))) == ["cs-calendar.csv", "cs-cells.csv", "cs-cohorts.csv", "cs-effects.csv",
            "cs-overall.csv", "cs-pretrends.csv", "panel-cohorts.csv", "panel-effects.csv", "panel-pretrends.csv"]
    end
    @testset "Callaway and Sant'Anna" begin
        options = (; unit = "unit", time = "year", outcome = "y", cohort = "first", distribution = "normal")
        estimate(table; changes...) = callaway_santanna(table; merge(options, changes)...)
        rounded(x) = round(x; digits = 6)
        pair(table, key, value) = (k = findfirst(==(value), getproperty(table, key)); rounded.((table.effect[k], table.se[k])))
        overall(r) = [pair(r.overall, :aggregation, a) for a in ("simple", "event time", "cohort", "calendar")]
        treated = map(c -> c[findall(!=(0), STAGGERED.first)], STAGGERED)

        # R's did 2.5.1, att_gt(bstrap = FALSE, cband = FALSE) and aggte(na.rm = TRUE): units never treated, varying base
        # period, doubly robust; each cohort's effect at a time, by event time, by cohort, by time and overall
        r = estimate(STAGGERED)
        cells = [(k = findfirst(j -> r.cells.cohort[j] == g && r.cells.time[j] == t, eachindex(r.cells.effect));
            rounded.((r.cells.effect[k], r.cells.se[k]))) for (g, t) in ((7, 5), (7, 8))]
        @test (cells, [pair(r.effects, :event_time, e) for e in (-2, 0, 3)], pair(r.cohorts, :cohort, 7), pair(r.cohorts, :cohort, 9),
            pair(r.calendar, :time, 7), pair(r.calendar, :time, 9), overall(r)) ==
            ([(0.010675, 0.089612), (1.996347, 0.170597)], [(-0.000241, 0.040001), (0.468011, 0.186927), (0.871297, 0.228423)],
            (1.987975, 0.175934), (-1.004267, 0.170014), (1.823025, 0.149978), (1.15345, 0.262189),
            [(1.378301, 0.178019), (1.843945, 0.112472), (1.004826, 0.227898), (1.293065, 0.135885)])
        # event times from -3 to 3 (min_e, max_e), and those of the cohorts seen 2 periods after treatment (balance_e)
        @test (overall(estimate(STAGGERED; leads = 3, lags = 3)), pair(estimate(STAGGERED; balance = 2).overall, :aggregation, "event time")) ==
            ([(0.670013, 0.20092), (0.670013, 0.20092), (0.670013, 0.20092), (1.293065, 0.135885)], (0.602918, 0.196206))
        # units not yet treated, universal base period, a covariate, weights and clusters coarser than the units: no
        # pre-trends test, as did, the covariance of the 14 cells before treatment being singular with 12 clusters
        r = estimate(STAGGERED; control = "notyet", base_period = "universal", covariates = ["x"], weights = "w", cluster = "region")
        @test ([pair(r.effects, :event_time, e) for e in (-3, 2)], overall(r), r.pretrends.df1, ismissing(r.pretrends.p)) ==
            ([(0.003828, 0.070708), (0.732151, 0.062883)],
            [(1.36834, 0.075432), (1.832825, 0.040267), (1.000777, 0.097618), (1.276919, 0.054015)], 14, true)
        # inverse probability weighting, and the outcome regression
        @test (pair(estimate(STAGGERED; method = "ipw", covariates = ["x"], weights = "w").effects, :event_time, 1),
            pair(estimate(STAGGERED; method = "ipw", covariates = ["x"], weights = "w").overall, :aggregation, "simple"),
            pair(estimate(STAGGERED; method = "reg", base_period = "universal", covariates = ["x"]).effects, :event_time, 2),
            pair(estimate(STAGGERED; method = "reg", base_period = "universal", covariates = ["x"]).overall, :aggregation, "simple")) ==
            ((0.591926, 0.199742), (1.367705, 0.165604), (0.736868, 0.195334), (1.371606, 0.160481))
        # no units never treated: the last cohort, before its treatment, and units not yet treated with a year's anticipation
        last = estimate(treated; covariates = ["x"])
        notyet = estimate(treated; control = "notyet", anticipation = 1)
        @test (pair(last.effects, :event_time, 0), pair(last.overall, :aggregation, "simple"), rounded(last.pretrends.statistic),
            pair(notyet.effects, :event_time, 0), pair(notyet.overall, :aggregation, "simple"), rounded(notyet.pretrends.statistic)) ==
            ((1.195519, 0.156497), (1.42411, 0.078378), 0.562784, (1.199719, 0.185252), (1.192297, 0.172061), 0.551189)
        @test_throws ArgumentError estimate(STAGGERED; method = "ols")

        # An unbalanced panel: half the units miss a year or two. R's did with allow_unbalanced_panel = TRUE, the cells
        # as repeated cross sections, doubly robust, then with units not yet treated, a universal base period, a
        # covariate, weights and clusters, inverse probability weighting and the outcome regression
        gaps = GAPS
        g7(r) = [(k = findfirst(j -> r.cells.cohort[j] == 7 && r.cells.time[j] == t, eachindex(r.cells.effect));
            rounded.((r.cells.effect[k], r.cells.se[k]))) for t in (5, 8)]
        r = estimate(gaps; unbalanced = true)
        @test (r.sample.panel, r.sample.unbalanced, r.sample.units, g7(r), [pair(r.effects, :event_time, e) for e in (-2, 0, 3)],
            pair(r.cohorts, :cohort, 7), overall(r), rounded.((r.pretrends.statistic, r.pretrends.p))) ==
            ("unbalanced", 30, 60, [(-0.036976, 0.098616), (1.992345, 0.175577)],
            [(0.009354, 0.046908), (0.494376, 0.195853), (0.923948, 0.235085)], (2.007894, 0.186063),
            [(1.404524, 0.181402), (1.869985, 0.119509), (1.031342, 0.230267), (1.319576, 0.139607)], (12.940849, 0.531186))
        r = estimate(gaps; unbalanced = true, control = "notyet", base_period = "universal", covariates = ["x"], weights = "w", cluster = "region")
        @test (g7(r), [pair(r.effects, :event_time, e) for e in (-2, 0, 3)], overall(r)) ==
            ([(-0.018657, 0.042311), (1.92571, 0.079773)], [(0.019522, 0.034192), (0.450209, 0.138436), (0.914768, 0.143759)],
            [(1.368598, 0.103154), (1.827711, 0.075404), (1.00824, 0.130123), (1.277964, 0.07789)])
        @test (overall(estimate(gaps; unbalanced = true, method = "ipw", covariates = ["x"], weights = "w")),
            overall(estimate(gaps; unbalanced = true, method = "reg", base_period = "universal", covariates = ["x"])),
            overall(estimate(gaps; unbalanced = true, covariates = ["x"], weights = "w"))) ==
            ([(1.386177, 0.182736), (1.852435, 0.11517), (1.016488, 0.237549), (1.294934, 0.136935)],
            [(1.396577, 0.175305), (1.860592, 0.106522), (1.025287, 0.224483), (1.314424, 0.134898)],
            [(1.374499, 0.164888), (1.831807, 0.091594), (1.013483, 0.216803), (1.284355, 0.121288)])
        # Without unbalanced the units with gaps are left out; with it, a balanced panel is as before
        @test ((s = estimate(gaps).sample; (s.panel, s.unbalanced, s.units)),
            (b = estimate(STAGGERED; unbalanced = true); (b.sample.panel, isequal(b.overall, estimate(STAGGERED).overall)))) ==
            (("balanced", 30, 30), ("balanced", true))

        # The multiplier bootstrap, against did's (bstrap = TRUE, cband = TRUE, biters = 100000): its standard errors
        # and simultaneous critical values, within the noise of the draws
        boot = estimate(STAGGERED; bootstrap = true, draws = 20000, seed = 1)
        c = boot.sample.critical
        @test all(isapprox.((c.cells, c.effects, c.cohorts, c.calendar), (2.797, 2.747, 2.319, 2.511); atol = 0.05))
        @test all(isapprox.((pair(boot.effects, :event_time, 0)[2], boot.overall.se[2], boot.overall.se[3]),
            (0.1886, 0.1138, 0.2297); rtol = 0.03))
        @test isapprox(boot.effects.high[1] - boot.effects.effect[1], c.effects * boot.effects.se[1]) && round(c.overall; digits = 6) == 1.959964
        @test boot.cells.effect == estimate(STAGGERED).cells.effect &&
            boot.effects.se == estimate(STAGGERED; bootstrap = true, draws = 20000, seed = 1).effects.se
        @test_throws ArgumentError estimate(STAGGERED; uniform = true)
    end

    @testset "stepped wedge" begin
        data(name) = CSV.File(joinpath(@__DIR__, "data", name))
        example, aggregated, growing = data("sw_data_example.csv"), data("sw_binomial.csv"), data("sw_growing.csv")
        design = (; cluster = "cluster", time = "period", treatment = "trt")
        estimate(table; options...) = (r = stepped_wedge(table; design..., options...); [r.estimate, r.se])

        # R's steppedwedge 1.0.0, analyze(method = "mixed"): its simulated sw_data_example (15 clusters, 5 sequences, 6
        # times), and a simulation of 18 clusters over 7 times whose effect grows with exposure time, where the immediate
        # treatment's effect is a third of the effects averaged over exposure times. The gaussian fits by REML
        @test isapprox(reduce(vcat, [estimate(example; outcome = "outcome_cont"),
                estimate(example; outcome = "outcome_cont", effect = "ETI"),
                estimate(example; outcome = "outcome_cont", effect = "ETI", estimand = "PTE", exposure = 3),
                estimate(growing; outcome = "outcome"),
                estimate(growing; outcome = "outcome", effect = "ETI"),
                estimate(growing; outcome = "outcome", effect = "ETI", estimand = "PTE", exposure = 4),
                estimate(growing; outcome = "outcome", effect = "ETI", exposure = 2:5)]),
            [0.0224086364381, 0.0741506876685, 0.0227260030266, 0.100722230828, -0.00517471111625, 0.117536641257,
                0.181421239158, 0.119627433432, 0.693272368616, 0.173459935563, 0.711445276615, 0.214323424958,
                0.753238974406, 0.171839087731]; atol = 1e-6)
        # The binomial fits by the Laplace approximation, of the 0 or 1 outcome and of the successes out of trials at
        # each cluster and time, which have the same likelihood (R's aggregated fit stops further from its maximum, so
        # both are compared with its fit of the 0 or 1 outcome), and the effect heterogeneous by exposure time. Where the
        # likelihood is flat, as with a variance near zero, R stops sooner, with a likelihood no higher than here; the
        # standard error of a PTE under TEH is that of the exposure time asked for (R's that of the one before)
        binomial = (; outcome = "successes", trials = "trials", family = "binomial")
        @test isapprox(reduce(vcat, [estimate(example; outcome = "outcome_bin", family = "binomial"),
                estimate(aggregated; binomial...),
                estimate(example; outcome = "outcome_bin", family = "binomial", effect = "ETI"),
                estimate(aggregated; binomial..., effect = "ETI"),
                estimate(example; outcome = "outcome_bin", family = "binomial", effect = "ETI", estimand = "PTE", exposure = 3),
                estimate(growing; outcome = "outcome", effect = "TEH"),
                estimate(growing; outcome = "outcome", effect = "TEH", estimand = "PTE", exposure = 4)]),
            [-0.189080063912, 0.132038210054, -0.189080063912, 0.132038210054, -0.183811820439, 0.159154613978,
                -0.183811820439, 0.159154613978, -0.304864523009, 0.19836704033, 0.524479642319, 0.185010235514,
                0.551594974196, 0.220107543027]; atol = 5e-5)
        @test_throws ArgumentError estimate(example; outcome = "outcome_cont", estimand = "PTE")
        @test_throws ArgumentError estimate(example; outcome = "outcome_cont", trials = "trt")

        # A study folder: the clusters are the units, and a TATE's exposure times [first, last]
        folder = mktempdir()
        cp(joinpath(@__DIR__, "data", "sw_growing.csv"), joinpath(folder, "panel.csv"))
        write(joinpath(folder, "study.toml"), """
            title = "Growing"
            source = "panel"
            unit = "cluster"
            time = "period"
            treatment = "trt"
            outcome = "outcome"
            estimators = ["stepped_wedge"]

            [stepped_wedge]
            effect = "ETI"
            exposure = [2, 5]
            """)
        output = mktemp() do path, io
            redirect_stdout(() -> run_study(folder), io)
            close(io)
            read(path, String)
        end
        @test (occursin("Stepped wedge: 2520 observations of 18 clusters over 7 times, 6 sequences; gaussian; exposure time indicators (ETI)", output),
            occursin("TATE (exposure times 2 to 5): 0.753, se 0.172, 95% interval 0.416 to 1.090", output),
            sort(readdir(joinpath(folder, "results")))) == (true, true, ["sw-effects.csv", "sw-estimate.csv"])
    end
end
