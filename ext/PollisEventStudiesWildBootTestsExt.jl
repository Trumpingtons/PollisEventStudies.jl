# The wild cluster bootstrap of panel_event_study, through David Roodman's WildBootTests.jl (MIT)
module PollisEventStudiesWildBootTestsExt

using PollisEventStudies, WildBootTests

function PollisEventStudies.wild_test(R, y, X, units, clusters, w; draws, rng, auxiliary, interval)
    test = wildboottest(R, zeros(size(R, 1)); resp = y, predexog = X, feid = units, clustid = clusters, obswt = w,
        reps = draws, rng, auxwttype = Symbol(auxiliary), getci = interval, getplot = false)
    (; p = WildBootTests.p(test), interval = interval ? vec(WildBootTests.ci(test)) : nothing, draws = WildBootTests.reps(test))
end

end
