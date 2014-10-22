## Caculate Trace of Hat matrix in LMM (only for PLSGeneral type in MixedModels)

function Λ0(m::LinearMixedModel)
	vv = SparseMatrixCSC{Float64,Int}[]
	for i in 1:length(m.λ)
		nl = length(m.facs[i].pool)
		ll = m.λ[i]
		p = size(ll,1)
		if p == 1
			isa(ll,MixedModels.PDScalF) || error("1×1 λ section should be a PDScalF type")
			push!(vv,ll.s .* speye(nl))
		else
			sfll = sparse(full(ll.ch.UL))
		push!(vv,blkdiag([sfll for _ in 1:nl]...))
		end
	end
	Triangular(blkdiag(vv...),:L,false)
end

function traceHat(m::LinearMixedModel)
	X=m.X.m
	XtX = Symmetric(X'X,:L)
	Zt=MixedModels.zt(m)
	Λt=Λ0(m).data

	Ztc = CholmodSparse(Λt*Zt,0)
	L = cholfact(Ztc,1.,true)
	
	perm=L.Perm .+ one(eltype(L.Perm))
	RX=cholfact(MixedModels.symcontents(XtX),:L)
	RZX = Λt*Zt*X

	copy!(RZX,solve(L, RZX[perm,:], CHOLMOD_L))

	LAPACK.potrf!('L',BLAS.syrk!('L','T',-1.,RZX, 1.,
		copy!(RX.UL, MixedModels.symcontents(XtX))))

	Cl= solve(L,solve(L,Ztc,CHOLMOD_P),CHOLMOD_L) 
	#Cl= solve(L,Ztc,CHOLMOD_L) 
	Cr=Ac_ldiv_B(tril(RX.UL), X' - full(Ac_mul_B(Cl,RZX))' )

	rho=sumabs2(Cl.nzval)+sumabs2(vec(Cr))
	rho
end