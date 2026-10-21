using PollisEventStudies, Test, Random, Dates, CSV, Tables

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
            println(io, "$(e.accession),1,\"Company $(e.ticker), Inc.\",$(e.ticker),8-K,2.02,2026-10-19,$(iso(e.accepted)),$(e.during),$next_open,$status,https://www.sec.gov/")
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
                println(io, "$(iso(t)),regular,$(e.ticker),$stock,$stock,$stock,$stock,$volume")
                println(io, "$(iso(t)),regular,$BENCHMARK,$level,$level,$level,$level,0")
            end
        end
    end
    folder
end

const EVENTS = [
    (accession = "E1", ticker = "AAA", accepted = DateTime(2026, 10, 19, 16, 0, 30), during = "true", next_open = nothing, jump = 0.02, beta = 1.2),
    # After the close on Friday: the jump comes at Monday's open
    (accession = "E2", ticker = "BBB", accepted = DateTime(2026, 10, 16, 21, 30, 10), during = "false", next_open = DateTime(2026, 10, 19, 14, 30), jump = -0.03, beta = 0.8),
    (accession = "E3", ticker = "CCC", accepted = DateTime(2026, 10, 20, 15, 0, 0), during = "true", next_open = nothing, jump = 0.01, beta = 1.0),
]

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
        @test default.tests.bmp_events == 3 && default.tests.abnormal_volume_events == 3
        # Five times the volume for the first 60 of the 120 minutes: log(5) / 2
        @test abs(default.tests.abnormal_volume - log(5) / 2) < 0.1
        # A period that does not follow the next open has no curve from it
        @test results[2].from_open === nothing

        files = write_results(results, joinpath(folder, "results"))
        @test basename.(files) == ["default-events.csv", "default-from-the-filing.csv", "default-from-the-next-open.csv", "short-run-events.csv", "short-run-from-the-filing.csv"]
        @test length(CSV.File(files[2])) == 30 + 120 + 1

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

    @testset "events still capturing" begin
        study = read_study(synthetic_study(mktempdir(); events = EVENTS[1:1], status = "capturing"))
        # Ten minutes after the filing, the CARs are known to the last complete minute and no further
        default = analyse(study; now = EVENTS[1].accepted + Minute(10))[1]
        known = findall(!ismissing, default.curves[1])
        @test last(known) == 30 + 1 + 8
        @test default.from_filing.events[end] == 0
        @test ismissing(default.tests.car)
    end

    @testset "no events" begin
        folder = mktempdir()
        write(joinpath(folder, "study.toml"), "title = \"Empty\"\nalert = \"x\"\n")
        write(joinpath(folder, "events.csv"), "accession,cik,name,ticker,form,items,filed,accepted,during_session,next_open,status,url\n")
        results = analyse(read_study(folder))
        @test length(results) == 1 && results[1].period.name == "Default" && results[1].from_filing === nothing
        @test occursin("0 events with prices", sprint(report, results))
    end
end
