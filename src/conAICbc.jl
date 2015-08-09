function V0inv(m::LinearMixedModel)

    Zt=MixedModels.zt(m)
    Lt=SelectMixedModels.Λ0(m)
    X=m.X.m
    Vinv = eye(size(m)[1]) - Zt'*inv(Zt*Zt' + inv(Lt*Lt'))*Zt
    
    return(Vinv - Vinv*X*inv(X'*Vinv*X)*X'*Vinv)
end 

function sumAY(m::LinearMixedModel)
     Vi = V0inv(m)
     return( (m.y' * Vi, Vi*m.y, m.y' * Vi *m.y)) 
end


function MakeBlkDeriv(m::LinearMixedModel)
    
    vv = SparseMatrixCSC{Float64,Int}[]
    for i in 1:length(m.λ)
    
        nl = length(m.facs[i].pool)
        ll = m.λ[i]
        #p = size(ll,1)
        p = size(ll)[1]
        if p == 1
            isa(ll,MixedModels.PDScalF) || error("1×1 λ section should be a PDScalF type")
            push!(vv,m.Ztblks[i]'*speye(nl)*m.Ztblks[i])
        else
            if isa(ll, MixedModels.PDDiagF)
                for k=1:p
                    sfll = sparse(zeros(p,p))
                    sfll[k,k]=1.0
                    push!(vv,m.Ztblks[i]'*blkdiag([sfll for _ in 1:nl]...)*m.Ztblks[i])
                end
            else  
                for k=1:p
                    for l=k:p  
                        sfll = sparse(zeros(p,p))
                        if k==l    
                            sfll[l,k]=1.0
                            push!(vv,m.Ztblks[i]'*blkdiag([sfll for _ in 1:nl]...)*m.Ztblks[i])
                        else
                            sfll[l,k]=1.0
                            sfll[k,l]=1.0
                            push!(vv,m.Ztblks[i]'*blkdiag([sfll for _ in 1:nl]...)*m.Ztblks[i])
                        end
                    end
                end
            end      
        end
    end
    #Triangular(blkdiag(vv...),:L,false)
    vv
end

function yAWAy(m::LinearMixedModel, W::Array{SparseMatrixCSC{Float64,Int64},1} )
    p=length(W)
    e = m.y-m.μ
    ewe=zeros(p)

    for i in 1:p
       ewe[i] =  (e' * W[i] *e)[1]
    end

    return(ewe)
end





function conAICbc(mm::LinearMixedModel)


    llk = condll(mm)

    oldf=mm.f
    newf=checkzerovar(mm)
    
    if oldf ==newf
        m=mm
    else
        retrms = filter(x->Meta.isexpr(x,:call) && x.args[1] == :|, ModelFrame(newf,mm.mf.df).terms.terms)
        length(retrms) > 0 || return(SelectMixedModels.IC( 0.0 , Inf, "VBc", [0.0,0.0], Inf ))
        m=fit(lmm(newf,mm.mf.df))
    end

    Zt=MixedModels.zt(m)
    Lt=SelectMixedModels.Λ0(m)
    X=m.X.m
    
    try 
        Vinv = eye(size(m)[1]) - Zt'*inv(Zt*Zt' + inv(Lt*Lt'))*Zt
    catch 
        Vinv = eye(size(m)[1])
    end       

    if countnz(Vinv) == size(m)[1]  
        return(SelectMixedModels.IC( 0.0 , Inf, "VBc", [0.0,0.0], Inf ))
    else
        Vinv = eye(size(m)[1]) - Zt'*inv(Zt*Zt' + inv(Lt*Lt'))*Zt
    end

    A= Vinv - Vinv*X*inv(X'*Vinv*X)*X'*Vinv
    rho =size(m)[1]-trace(A)
     
   

    (nn,pp,qq,tt)=size(m)

    W=MakeBlkDeriv(m)
    eWe=yAWAy(m, W)
    e=m.y-m.μ
    tye= (m.y'*e)[1]
    theta=MixedModels.θ(m)
    s=length(theta)

    B=zeros(s,s)

    for j in 1:s 
        for l in 1:j 
            B[j,l] =  -tye^2 *trace(W[j]*Vinv*W[l]*Vinv)/nn 
            B[j,l] -=  eWe[j] * eWe[l] +2.0* (e' * W[l] *A* W[j] *e)[1]*tye 
            
            if j != l 
                B[j,l] =B[l,j]
            end
         end
    end

    Binv=inv(B)

    G= zeros(s,nn)

    
    for j in 1:s 
      G[j,:] = 2.0* (tye * (e' * W[j] *A - eWe[j] * e')[1] )
    end

     

    bc=0.0

    for j in 1:s 
      bc += sum ( Binv[j,:]*G*A*W[j]*e )
    end

    df = rho+bc+1.0
    val = -2.0*llk+2.0*df

    cr=[1.0,1.0]

    SelectMixedModels.IC(llk, df, "VBc", cr, val )
end

 