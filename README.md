# PollisEventStudies.jl

Intraday event studies of SEC filings, on the study folders [Pollis](https://pollis.org) keeps in `~/.pollis/event-studies`.

Pollis captures the minute prices of each filer's stock and of a benchmark around the filings an alert fires. Each study folder holds:

- `study.toml`: the benchmark, the model of normal returns (market-adjusted or the market model) and the periods, in minutes before and after the filing;
- `events.csv`: the filings, with the time EDGAR accepted each (UTC), whether it came in during the regular session and the next open after it;
- `prices/`: a file per filing with the minute bars of the stock and of the benchmark.

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
- at the end of the period: the CAR's cross-sectional t, the Boehmer, Musumeci and Poulsen (1991) t on standardized CARs, and abnormal volume (regular-session minutes after the filing against those before it).

Events still capturing prices count up to the last complete minute. The standardized CARs need the estimation days, which Pollis captures for the market model only.

## License

MIT
