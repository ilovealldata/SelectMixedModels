 ## Return a block in the Zt matrix from one term.
function ztblk(m::Matrix,v)
    nr,nc = size(m)
    nblk = maximum(v)
    NR = nr*nblk                        # number of rows in result
    cf = length(m) < typemax(Int32) ? int32 : int64 # conversion function
    SparseMatrixCSC(NR,nc,
                    cf(cumsum(vcat(1,fill(nr,(nc,))))), # colptr
                    cf(vec(reshape([1:NR],(nr,int(nblk)))[:,v])), # rowval
                    vec(m))            # nzval
end
ztblk(m::Matrix,v::PooledDataVector) = ztblk(m,v.refs)

##  Conditional Likelihood
function condll(m::LinearMixedModel)
    n,p,q,k = size(m)
    sig = std(m)[end][1]
    cll=0.0
    for i in 1:n
        cll += logpdf(Normal(m.μ[i],sig),m.y[i])
    end
    cll
end


function checkzerovar(m::LinearMixedModel)
    ff = deepcopy(m.f)
    # set fixed term=true, reterm=false in rhs of formula
    fixterms=map(x-> !(Meta.isexpr(x,:call) && x.args[1] == :|), m.mf.terms.terms)
    #find λ and their lowe bounds , compare them, 
    # set true if boundary conditions are met or set false otherwise                
    est=map(MixedModels.θ,m.λ)
    lowbd=map(MixedModels.lower,m.λ)
    remainterms1=map(all,Any[l .> k  for (l,k) in zip(est,lowbd)])

    # set reterm be false (will be excluded) if  boundary conditions are not met 
    fixterms[fixterms .== false] = remainterms1
    # chaneg rhs of formula by excluding reterms for which boundary conditions are not met 
    
    ffm=m.mf.terms.terms[convert(Array{Bool,1},fixterms)]
    
    if all(remainterms1)
        ff=ff
    elseif length(ffm)==0
        warn("Est. of variance components of all terms are on the boundary")
        ff=y~1
    elseif (length(ffm)==1 && !any(remainterms1))
        warn("Est. of variance components of all terms are on the boundary")    
        ff=eval(y~$(ffm[1]))
    elseif (length(ffm)==1 && any(remainterms1))    
        ff=eval(y~$(ffm[1]))
    # Note that 
    # if there are more than 1 term, the first element of formula.rhs.args is :+ with 
    # head function :call
    else    
        a=Expr(:call,tuple ( [:+,ffm]...)...)
        ff=eval(y~$a)
    end 

    #returm modofied formula
    ff
end