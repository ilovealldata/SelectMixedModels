type conAIC <: AIC
    condlike :: Float64 
    tracehat :: Float64
    method   :: ASCIIString
    corterm  :: Vector
    value    :: Float64
end


function conAIC(m::LinearMixedModel)

    rho = traceHat(m)
    llk = condll(m)

    (nn,pp,qq,tt)=size(m)
    
    cr=zeros(Float64,2)
    cr[1]=nn*(nn-pp-1.0)/((nn-pp)*(nn-pp-2.0))
    cr[2]=nn*(pp+1.0)/((nn-pp)*(nn-pp-2.0))
    
    conAIC(llk, rho, "VB", cr, -2.0*llk+2.0*(cr[1]*(1+rho)+cr[2]) )
end