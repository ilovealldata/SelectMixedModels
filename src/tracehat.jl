## Caculate Trace of Hat matrix in LMM (only for PLSGeneral type in MixedModels)

function traceHat(m::LinearMixedModel)

    typeof(m.s) ==PLSGeneral{Int32} || error("Method available only for PLSGeneral type")
	
	Cl= solve(m.s.L,solve(m.s.L,m.s.λtZt,CHOLMOD_P),CHOLMOD_L) 
	Cr=Ac_ldiv_B(tril(m.s.RX.UL), m.X.m' - full(Ac_mul_B(Cl,m.s.RZX))' )
	#(sumabs2(Cl.nzval)+sumabs2(vec(Cr)),Cl,Cr)
	rho=sumabs2(Cl.nzval)+sumabs2(vec(Cr))
    rho

end