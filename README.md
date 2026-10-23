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

- the average abnormal return (AAR) and the cumulative average abnormal return (CAAR) a minute at a time, from the filing and, for filings outside the regular session, from the next open, each with its cross-sectional t and the CAAR's 95% band (`caar_lower`, `caar_upper`, where there are two events or more), and each event's CAR a minute at a time (`curves`, and `open_curves` from the next open);
- at the end of the period: the CAR's cross-sectional t, the Boehmer, Musumeci and Poulsen (1991) t on standardized CARs, the same t adjusted for the correlation of the events' abnormal returns (Kolari and Pynn&ouml;nen 2010), the generalized rank test (Kolari and Pynn&ouml;nen 2011), Cowan's (1992) generalized sign test, the Wilcoxon signed-rank test, abnormal volume (regular-session minutes after the filing against those before it) and abnormal volatility (the log of the realized volatility of the minute abnormal returns after the filing over that of as many regular-session minutes before it), each with its two-sided p-value.

A filing outside the regular session, in a period that follows the next open, is tested at the end of the period's minutes after the next open, where the market first trades on it. Events still capturing prices count up to the last complete minute. The standardized CARs need the estimation days, which Pollis captures for the market model only.

The adjustment for cross-correlation multiplies the BMP t by sqrt((1 - r) / (1 + (n - 1) r)), r the mean correlation over all pairs of the n events: the correlation of a pair's standardized abnormal returns over the estimation minutes they share when their periods overlap, and zero when they do not. The generalized rank test ranks each event's standardized CAR, scaled by the CARs' cross-sectional standard deviation, among its standardized abnormal returns on the estimation days. The generalized sign test sets how many CARs are positive against how many the share of positive abnormal returns on the estimation days leads one to expect; the Wilcoxon signed-rank test, by [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl), sums the ranks of the positive CARs among their absolute values. Both rest on signs and ranks rather than means, for heavy-tailed returns and few events. The p-values of the t tests come from Student's t distribution, that of the sign test from the normal, through [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

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
- the mean move a minute at a time from the release, with its 95% band (`paths`, by instrument), and each release's (`curves`, by instrument).

`write_results` writes them to `<period>-moves.csv`, `<period>-summary.csv` and `<period>-paths.csv`.

## Panel event studies

Units of a panel first treated at different times, such as firms adopting a rule, states passing a law or hospitals
joining a stepped-wedge trial, against units not yet or never treated:

```julia
results = panel_event_study(table; unit = "state", time = "year", outcome = "y", cohort = "first_treated",
    leads = 5, lags = 5)
report(results)
```

`table` is any Tables.jl table, a row a unit and a time, with times in whole numbers (years, or quarters or months
numbered). When each unit was first treated comes from `cohort`, a column of that time (missing or 0 for never), or from
`treatment`, the first time its column is not zero. Event times run from `-leads` to `lags` against `reference` (-1);
those beyond the ends count at the ends. Units treated from the first time are left out.

- **Two-way fixed effects** (TWFE): the outcome on unit and time effects and an indicator of each event time. With
  effects that differ from cohort to cohort, its effects mix in those of other event times, before treatment too.
- **Sun and Abraham (2021)**: each cohort's effect at each event time against the control units (`control = "never"`,
  the units never treated, or `"last"`, the units treated last, before their treatment), averaged over the cohorts
  with their shares of the event time's observations.

Standard errors are clustered by `cluster` (the unit by default), with Stata's small-sample factor unless
`small_sample = false`. p-values and the 95% bands come from t with the clusters less one degrees of freedom
(`distribution = "t"`, as Stata and fixest) or from the normal (`distribution = "normal"`, as did). Sun and Abraham's variance adds that of the
cohorts' shares, clustered, unless `fixed_shares = true` holds them fixed. Each result has the effects by event time (`effects`), each cohort's (`cohorts`), the Wald
tests that the effects before treatment are zero (F under t, chi-squared under the normal) (`pretrends`) and the sample. `weights` weighs the observations, and
`covariates` are added to both regressions.

To compare with R:

| R | Keywords |
|---|---|
| fixest, `sunab` and `i()` | `fixed_shares = true` |
| did, `att_gt` (units never treated, universal base period) and `aggte(type = "dynamic")` | `small_sample = false, distribution = "normal"` (did's `bstrap = FALSE, cband = FALSE`) |

On a test panel with weights, a covariate and clusters coarser than the units, the cohorts' effects, the TWFE
effects, the effects by event time and their standard errors are those of fixest to 6 decimals; on one without, the
effects by event time and their standard errors, the variance of the shares included, their p-values and their 95% bands are those of did to 6 decimals. Stata's eventstudyinteract takes the
shares' variance observation by observation.

With `uniform = true`, the 95% bands of each estimator's effects are simultaneous: they cover all the event times
together 95% of the time, as an event-study plot needs. The critical value is the sup-t of Montiel Olea and
Plagborg-Moller (2019): the 95% quantile of the largest |t| over the event times, drawn (`draws`, 10,000 by default;
`seed` to repeat them) from the effects' correlations, normal or Student's t as `distribution`. The report gives each
estimator's. On the test panel, 2.777 for Sun and Abraham's 16 event times under the normal, between the pointwise
1.960 and Bonferroni's 2.955. With `multiplier = true`, the critical value comes from did's multiplier bootstrap
instead, on each cluster's influence on the effects, each weighed by +1 or -1 in each draw: 2.706 with 200,000 draws,
against did's 2.708. Over 60 clusters, weights of +1 or -1 give the largest |t| thinner tails than normal draws do.

With few clusters, clustered standard errors are unreliable. With `wild = true`, after `using WildBootTests`, each
effect also has the p-value of the wild cluster bootstrap (Cameron, Gelbach and Miller 2008; Roodman, MacKinnon, Nielsen
and Webb 2019), with the null imposed, and its 95% interval by inverting the test, and the pre-trends test has its own
p-value. Each draw weighs each cluster's residuals by +1 or -1 (`wild_weights = "rademacher"`, every pattern of signs
when there are no more than `draws`) or by Webb's six values (`"webb"`, better with fewer than about 12 clusters). Sun and
Abraham's effects hold the cohorts' shares fixed. The bootstrap runs in David Roodman's
[WildBootTests.jl](https://github.com/droodman/WildBootTests.jl) (MIT), loaded only when asked for. Against R's
fwildclusterboot (`boottest` on fixest's models), with 12 clusters and so all 4,096 patterns of signs, the p-values are
equal to 6 decimals; the intervals' bounds differ by up to 0.0005, both programs stopping on the step of the p-value
past the bound (TWFE's upper bound at event time 3: 1.158875, R 1.158985, WildBootTests 1.159517).

A study folder with `source = "panel"` in `study.toml` names the table (`data`, `panel.csv` by default), the same
keywords and the estimators, `panel_event_study` (the default), `callaway_santanna` and `stepped_wedge` (below), any of
them, whose own keywords go in a `[callaway_santanna]` or a `[stepped_wedge]` table:

```toml
title = "Minimum wage"
source = "panel"
data = "panel.csv"
unit = "state"
time = "year"
outcome = "y"
cohort = "first_treated"
leads = 5
lags = 5
estimators = ["panel_event_study", "callaway_santanna"]

[callaway_santanna]
control = "notyet"
method = "dr"
```

`run_study` reports each and writes `panel-effects.csv`, `panel-cohorts.csv` and `panel-pretrends.csv`, and the
`cs-*.csv` files of Callaway and Sant'Anna.

### Callaway and Sant'Anna (2021)

Each cohort's effect at each time, ATT(g, t), a difference in differences against the control units, and its
aggregations:

```julia
results = callaway_santanna(table; unit = "state", time = "year", outcome = "y", cohort = "first_treated",
    control = "notyet", covariates = ["x"])
report(results)
write_results(results, "results")   # cs-cells.csv, cs-effects.csv, cs-cohorts.csv, cs-calendar.csv, cs-overall.csv, cs-pretrends.csv
```

- `control`: `"never"`, `"notyet"` (never treated and not yet treated), `"last"` (never treated and the cohort treated
  last, before its treatment) or `"auto"` (`"never"` when there are units never treated, else `"last"`).
- `base_period`: `"varying"` (the time before, for the cells before treatment) or `"universal"` (the last time before
  the cohort's treatment for all).
- `method`: `"dr"`, doubly robust (Sant'Anna and Zhao 2020), `"ipw"` or `"reg"`, with `covariates` taken at the earlier
  of the two times.
- `anticipation`, the times before its treatment a cohort may respond; `weights` and `cluster` as above.

The cells are aggregated by event time (from `-leads` to `lags`; `balance = e` keeps the cohorts seen e times after
treatment), by cohort, by calendar time and overall, weighed by the cohorts' sizes, the standard errors counting the
estimation of those weights. Standard errors come from the influence functions, clustered, and the test of pre-trends
is the Wald test that the cells before treatment are zero (none when their covariance is singular, as in did).

p-values and 95% bands come from the normal, as did, the authors' own package (`distribution = "normal"`, the
default), or from Student's t with the clusters less one degrees of freedom (`distribution = "t"`), wider with few
clusters.

The code is a port of R's [did](https://github.com/bcallaway11/did) (`att_gt`, `aggte`) and of
[DRDID](https://github.com/pedrohcgs/DRDID) (`drdid_panel`, `std_ipw_did_panel`, `reg_did_panel`). With the normal
(the default), against did's `bstrap = FALSE, cband = FALSE`, on the test panel, with and without units
never treated: the cells, the effects by event time, by cohort and by calendar time, the overall averages, their
standard errors and the test of pre-trends are those of did to 6 decimals, for each control group, base period and
method, with weights, a covariate, coarser clusters, anticipation, `min_e`/`max_e` (`leads`/`lags`) and `balance_e`
(`balance`).

With `bootstrap = true`, the standard errors come from did's multiplier bootstrap (`draws`, 1000 by default, each
weighing each cluster's influence by +1 or -1; `seed` to repeat them), and the 95% bands of the cells and of the effects
by event time, cohort and time are simultaneous (`uniform`, did's `cband`): they cover all of a table's rows together
95% of the time, as an event-study plot needs. The report gives each table's critical value. Against did with 100,000
draws, the critical values and standard errors agree within the noise of the draws (2.743 and 2.747 for the event times).
Units observed at only some times are left out, as in did, unless `unbalanced = true` (did's
`allow_unbalanced_panel`): each ATT(g, t) then compares the observations at t and at the base period as repeated cross
sections, with DRDID's `drdid_rc`, `std_ipw_did_rc` and `reg_did_rc`, each observation with its own covariates and
weight, and a unit's influence is the sum of its observations'. On a panel where half the units miss a year or two, it
agrees with did to 6 decimals in five configurations (each method, the not-yet-treated control, a universal base period,
weights, a covariate and coarser clusters). `panel_event_study` needs no option: its regression takes every observation,
and it agrees with fixest on the same panel.

### Stepped wedges

When the times clusters start a treatment are drawn at random, a stepped-wedge cluster randomized trial, the mixed
model of Hussey and Hughes (2007): a fixed effect for each time, a random intercept for each cluster and one for each
cluster at each time. After `using MixedModels`:

```julia
using MixedModels
result = stepped_wedge(table; cluster = "clinic", time = "period", treatment = "trt", outcome = "y", effect = "ETI")
result.estimate, result.se, result.interval   # the effect averaged over exposure times (TATE)
result.by_exposure                              # the effect at each exposure time
```

- `effect`: `"IT"`, one effect at every exposure time (the time since the cluster's first treatment, plus 1); `"ETI"`,
  an effect at each exposure time; `"TEH"`, a random effect at each exposure time about a mean effect. When the effect
  grows with exposure time, IT's estimate can be far from the average effect (Kenny et al. 2022): in the test
  simulation, 0.18 against ETI's 0.69.
- `estimand`: `"TATE"`, the effects averaged over the exposure times `exposure` (a range, all by default), or `"PTE"`,
  the effect at exposure time `exposure`.
- `family`: `"gaussian"` (REML, as lme4's `lmer`) or `"binomial"` (Laplace, as `glmer`), of a 0 or 1 outcome or of
  successes out of `trials`; `individual` adds a random intercept for each person.

The result also has the standard deviations of the random effects, the intracluster correlation of a gaussian model and
the fitted MixedModels.jl model; `report` prints it and `write_results` writes `sw-effects.csv` (by exposure time) and
`sw-estimate.csv`. `PollisEventStudies.result_tables(result)` gives the tables `write_results` writes for any of the
panel estimators' results, each with its file's name, without writing them. In a study folder, `estimators = ["stepped_wedge"]`; the clusters are `cluster`, or else `unit`, and
a TATE's exposure times go in the `[stepped_wedge]` table as `exposure = [2, 5]`:

```toml
source = "panel"
unit = "clinic"
time = "period"
treatment = "trt"
outcome = "y"
estimators = ["stepped_wedge"]

[stepped_wedge]
effect = "ETI"
exposure = [2, 5]
```

The code is a port of the mixed-model part of R's
[steppedwedge](https://cran.r-project.org/package=steppedwedge) 1.0.0 (`analyze(method = "mixed")`, calendar time
categorical), fitted by [MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl) (MIT), loaded only when asked
for. On steppedwedge's simulated example and on a simulation whose effect grows with exposure time, the gaussian IT and
ETI estimates and standard errors are those of R within 1e-6. The binomial fits and TEH are within 5e-5: where a
variance is near zero the likelihood is flat and R stops sooner, at a likelihood no higher than here
(-1376.65973195 against -1376.65973191). Under TEH, the standard error of a PTE takes the conditional variance of the
exposure time asked for; steppedwedge takes that of the exposure time before it (0.220108 against R's 0.213374 at
exposure time 4).

## Vignette

[Event Studies: a filing's minutes, a law's years](vignette/PollisEventStudies.pdf)

## License

MIT
