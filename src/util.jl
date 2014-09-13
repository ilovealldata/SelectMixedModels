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


