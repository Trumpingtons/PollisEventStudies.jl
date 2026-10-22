# PollisEventStudies.jl

Intraday event studies of SEC filings and of the FOMC's releases, on the study folders [Pollis](https://pollis.org) keeps in `~/.pollis/event-studies`.

Pollis captures the minute prices of each filer's stock and of a benchmark around the filings an alert fires. Each study folder holds:

- `study.toml`: the benchmark, the model of normal returns (market-adjusted or the market model) and the periods, in minutes before and after the filing;
- `events.csv`: the filings, with the time EDGAR accepted each (UTC) and the event's time (`event_time`: the earlier of the newswire release and the acceptance, or the acceptance alone, as `study.toml` says; the acceptance when the column is not there), whether the event came during the regular session and the next open after it;
- `prices/`: a file per filing with the minute bars of the stock and of the benchmark.

For the FOMC's releases (`source = "fomc"` in `study.toml`), Pollis captures the bars of the study's instruments around each release in the Federal Reserve's Monetary Policy Press Releases feed: federal funds futures (FF1 to FF4), SOFR futures (SR1, SFR1 to SFR5), the 2-year note future (ZT), the 5- and 10-year yields (FVX, TNX), TIPS and the S&P 500. The folder holds:

- `study.toml`: the releases (statement, projections, minutes, other), the instruments and the periods, in minutes before and after the release;
- `events.csv`: the releases, with their time (UTC) and the minutes each bar spans (1, or 5 for a release older than 30 days);
- `instruments.csv`: for each release, each instrument's Yahoo Finance symbol (the contract months follow from the release's month), its unit and its scale;
- `prices/`: a file per release with the bars of each symbol.

## Usage

```julia
using PollisEventStudies

results = run_study("path/to/study")   # prints the report and writes results/ in the folder
```

or step by step:

```julia
study = read_study("path/to/study")
results = analyse(study)               # one result per period
report(results)
write_results(results, "results")
```

For each period:

- the average abnormal return (AAR) and the cumulative average abnormal return (CAAR) a minute at a time, from the filing and, for filings outside the regular session, from the next open, each with its cross-sectional t;
- at the end of the period: the CAR's cross-sectional t, the Boehmer, Musumeci and Poulsen (1991) t on standardized CARs, the same t adjusted for the correlation of the events' abnormal returns (Kolari and Pynn&ouml;nen 2010), the generalized rank test (Kolari and Pynn&ouml;nen 2011), and abnormal volume (regular-session minutes after the filing against those before it), each with its two-sided p-value.

A filing outside the regular session, in a period that follows the next open, is tested at the end of the period's minutes after the next open, where the market first trades on it. Events still capturing prices count up to the last complete minute. The standardized CARs need the estimation days, which Pollis captures for the market model only.

The adjustment for cross-correlation multiplies the BMP t by sqrt((1 - r) / (1 + (n - 1) r)), r the mean correlation over all pairs of the n events: the correlation of a pair's standardized abnormal returns over the estimation minutes they share when their periods overlap, and zero when they do not. The generalized rank test ranks each event's standardized CAR, scaled by the CARs' cross-sectional standard deviation, among its standardized abnormal returns on the estimation days. The p-values come from Student's t distribution, through [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

The same for groups of the events, with `group_by = "form"`, `"item"` (each 8-K item; a filing with several is in the group of each), `"session"` (during the regular session or outside it), `"ticker"`, or a function of the event that gives a group's name or several, in `study.toml` or as a keyword:

```julia
results = analyse(study; group_by = "item")
results[1].groups                      # each with its name, events, AAR and CAAR, and tests
run_study("path/to/study"; group_by = event -> event.form == "8-K" ? "8-K" : "Other")
```

`write_results` then also writes `<period>-groups.csv` (a group a row) and the groups' AAR and CAAR, each row with its group.

A stock that does not trade every minute has no bar in the minutes it did not trade in. The market model's returns are measured from one trade to the next, against the benchmark's move over the same minutes, and weighed by the minutes they span; a return of zero in each minute without a trade would bias beta towards zero. The standardized CAR's variance counts every minute of the regular session in the period, whether the stock traded in it or not.

Grouping the events and the treatment of thin trading follow [EventStudy.jl](https://github.com/florian-st/EventStudy.jl) (MIT), whose thin-trading treatment is that of Stata's eventstudy2.

For each period of a study of FOMC releases:

- each instrument's move from the start of the period to its end at each release (`moves`, a release a row): in basis points of the implied rate (100 less the price) for the federal funds and SOFR futures, of the yield for the yields, and in per cent (log) for prices; FF1's move is scaled by D / (D - d), D the days of the release's month and d its day (Kuttner 2001), and in the month's last week FF2's contract stands in for it, unscaled;
- the mean move of each instrument with its cross-sectional t (`summary`);
- the mean move a minute at a time from the release (`paths`, by instrument).

`write_results` writes them to `<period>-moves.csv`, `<period>-summary.csv` and `<period>-paths.csv`.

## License

MIT
