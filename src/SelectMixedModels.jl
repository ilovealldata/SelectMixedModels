module SelectMixedModels

	using MixedModels, ArrayViews, DataArrays, DataFrames, Distributions, PDMats, StatsBase 
    using Base.LinAlg.CHOLMOD: CholmodFactor, CholmodSparse, CholmodSparse!,
        chm_scale, CHOLMOD_SYM, CHOLMOD_L, CHOLMOD_Lt, CHOLMOD_P, CHOLMOD_ROW, solve
	using Base.LinAlg: Cholesky, Ac_ldiv_B!, A_rdiv_Bc!, chksquare, transpose!


export
    conAIC,		# Conditional AIC
    marAIC,     # marginal AIC
    marBIC,     # marginal BIC
    lmmg,       # fit a linear mixed-effects model (LMM) 
                # with PLSGeneral type only
    condll,     # compute conditional likelihood (LMM)
    traceHat    # Caculate Trace of Hat matrix in LMM

    include("lmmgeneral.jl")
    include("tracehat.jl")
    include("lmmgeneral.jl")
    include("conAIC.jl")
    include("marIC.jl")
    include("util.jl")

end # module
