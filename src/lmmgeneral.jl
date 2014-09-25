 function lmmg(f::Formula, fr::AbstractDataFrame)
    mf = ModelFrame(f,fr)
    X = ModelMatrix(mf)
    y = convert(Vector{Float64},DataFrames.model_response(mf))
    Xty = X.m'y

    retrms = filter(x->Meta.isexpr(x,:call) && x.args[1] == :|, mf.terms.terms)
    length(retrms) > 0 || error("Formula $f has no random-effects terms")

    grps = {t.args[3] for t in retrms}       # expressions for grouping factors
    facs = {pool(getindex(mf.df,grp)) for grp in grps}
    l = Int[length(f.pool) for f in facs]
    if (perm = sortperm(l;rev=true)) != [1:length(grps)]
        permute!(retrms,perm)
        permute!(grps,perm)
        permute!(facs,perm)
        permute!(l,perm)
    end
    Xs = {MixedModels.lhs2mat(t,mf.df)' for t in retrms} # transposed model matrices
    p = Int[size(x,1) for x in Xs]
    λ = {pp == 1 ? MixedModels.PDScalF(1.,1) : 
                MixedModels.PDLCholF(cholfact(eye(pp),:L)) for pp in p}
    if length(unique(grps)) < length(grps)
        grps,Xs,p,λ,facs,l = MixedModels.amalgamate(grps,Xs,p,λ,facs,l)
    end
    q = sum(p .* l)
    uβ = zeros(q + size(X,2))
    Zty = {zeros(pp,ll) for (pp,ll) in zip(p,l)}
    u = {}
    offset = 0
    for (x,ff,zty) in zip(Xs,facs,Zty)
        push!(u,contiguous_view(uβ, offset, size(zty)))
        offset += length(zty)
        for (j,jj) in enumerate(ff.refs)
            for i in 1:size(zty,1)
                zty[i,jj] += y[j] * x[i,j]
            end
        end
    end
    local s
    
    Zt = vcat(map(ztblk,Xs,facs)...)
    s =  PLSGeneral(Zt,X.m,facs)
    
    LinearMixedModel(false, X, Xs, Xty, map(ztblk,Xs,facs), Zty,
                     map(zeros, u), f, facs, false, map(string,grps),
                    [zeros(k,k) for k in p], mf, similar(y),
                     s, u, uβ, y, λ, similar(y))
end



