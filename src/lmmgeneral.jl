 function lmmg(f::Formula, fr::AbstractDataFrame)
    mf = ModelFrame(f,fr)
    X = ModelMatrix(mf)
    y = convert(Vector{Float64},DataFrames.model_response(mf))
    Xty = X.m'y
                                        # process the random-effects terms
    retrms = filter(x->Meta.isexpr(x,:call) && x.args[1] == :|, mf.terms.terms)
    length(retrms) > 0 || error("Formula $f has no random-effects terms")
    Xs = [lhs2mat(t,mf.df)' for t in retrms] # transposed model matrices
    p = Int[size(x,1) for x in Xs]
    λ = Any[pp == 1 ? PDScalF(1.,1) : PDLCholF(cholfact(eye(pp),:L)) for pp in p]
    grps = Any[t.args[3] for t in retrms]
    ugrps = unique(grps)
    if length(ugrps) < length(grps)     # amalgamate terms with the same grouping factor
        for g in ugrps
            ii = [1:length(grps)][grps .== g]
            length(ii) == 1 && continue
            iii = copy(ii)
            i1 = shift!(ii)
            p[i1] = sum(p[iii])
            deleteat!(p,ii)
            Xs[i1] = vcat(Xs[iii]...)
            deleteat!(Xs,ii)
            λ[i1] = amalgamate(λ[iii])
            deleteat!(λ,ii)
            deleteat!(grps,ii)
        end
    end
    facs = [getindex(mf.df,g) for g in ugrps]
    map!(x->isa(x,PooledDataArray) ? x : pool(x), facs)
    l = Int[length(f.pool) for f in facs]

    q = sum(p .* l)
    uβ = zeros(q + size(X,2))
    Zty = [zeros(pp,ll) for (pp,ll) in zip(p,l)]
    u = Any[]
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



