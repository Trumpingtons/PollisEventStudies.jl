# The stepped-wedge mixed model of stepped_wedge, through MixedModels.jl (MIT)
module PollisEventStudiesMixedModelsExt

using PollisEventStudies, MixedModels
using MixedModels.StatsModels: term

function PollisEventStudies.mixed_fit(data, fixed, groups, slope; family, weights)
    random = [term(1) | term(g) for g in groups]
    slope === nothing || push!(random, (term(0) + term(slope)) | term(:exposure))
    formula = term(:y) ~ foldl(+, [term(0); term.(fixed); random])
    contrasts = Dict{Symbol, Any}(g => Grouping() for g in [groups; :exposure])
    model = if family == "gaussian"
        LinearMixedModel(formula, data; contrasts)
    elseif weights === nothing
        GeneralizedLinearMixedModel(formula, data, Bernoulli(); contrasts)
    else
        GeneralizedLinearMixedModel(formula, data, Binomial(); contrasts, wts = weights)
    end
    # The default tolerances stop short where a variance is near zero, as that of the clusters at a time often is: the
    # tighter ones reach lme4's likelihood
    model.optsum.ftol_rel = 1e-14
    model.optsum.ftol_abs = 1e-12
    family == "gaussian" ? fit!(model; REML = true, progress = false) : fit!(model; progress = false)
    sd = Dict{Symbol, Float64}(g => only(values(s)) for (g, s) in pairs(model.sigmas))
    modes = Dict{String, Tuple{Float64, Float64}}()
    if slope !== nothing
        k = findfirst(re -> MixedModels.fname(re) == :exposure, model.reterms)
        re = model.reterms[k]
        for (level, mode, variance) in zip(re.levels, vec(ranef(model)[k]), vec(condVar(model)[k]))
            modes[string(level)] = (mode, variance)
        end
    end
    (; names = Symbol.(coefnames(model)), coef = coef(model), vcov = vcov(model), sd, modes,
        residual = family == "gaussian" ? model.sigma : nothing, model)
end

end
