function conAIC(m::LMMScalar1)
  q=length(m.L); theta=m.theta; p,n = size(m.Xt);
  ## note
  ## L*Rzx =lambdat *ZtX  => Rzx =theta*L^-1 *ZtX
  ## Hl = [ L ,0 ; Rzx , Rx]
  Hlinv=inv(tril(hvcat((2,2),diagm(m.L),(theta ./ m.L).*m.XtZ',(theta ./ m.L)'.*m.XtZ, m.RX.UL')))
  H= Hlinv' * Hlinv  
  #making Zt
  Zt=zeros(Float64,q,n)
  for i in 1:n 
    j = m.Ztrv[i]; z = m.Ztnz[i]
    Zt[j,i] = z 
  end
  # rho= tr ( [Z*lambda; X] H [lambda] [lambdat*Zt, Xt] )
  rho = trace( hcat(theta*Zt',m.Xt') * H * vcat(theta*Zt, m.Xt))
  #REML
  #rho=(n-p-1.0)*(rho+1.0)/(n-p-2.0) + (p+1.0)/(n-p-2.0)
  #MLE
  rho=n/(n-p-2.0) * ((rho+1) -(rho-p)/(n-p))
  (objectivecon(m)+2.0*rho, rho) 
end

function conAICb(m::LinearMixedModel)
  n,p,q,k = size(m)
  spλtZt=sparse(m.s.λtZt)
  k12= spλtZt * sparse(m.X.m)  
  #Hlinv = inv(hvcat((2,2),m.s.λtZt *m.s.λtZt' +eye(q) ,k12,k12',sparse(m.X.m)'*sparse(m.X.m)))
  Hlinv = inv(tril(full([ spλtZt *spλtZt'+speye(q) k12 ; k12' sparse(m.X.m)'*sparse(m.X.m)])))
  rho = trace( vcat(spλtZt,sparse(m.X.m)')*vcat(spλtZt,sparse(m.X.m)')' *Hlinv)
  #REML
  #rho=(n-p-1.0)*(rho+1.0)/(n-p-2.0) + (p+1.0)/(n-p-2.0)
  #MLE
  rho=n/(n-p-2.0) * ((rho+1) -(rho-p)/(n-p))
  (-2.0*condll(m)+2.0*rho, rho) 
end

function conAICa(m::LinearMixedModel)
  n,p,q,k = size(m)
  #k12= m.s.λtZt * sparse(m.X.m)  
  #Hlinv = inv(hvcat((2,2),m.s.λtZt *m.s.λtZt' +eye(q) ,k12,k12',sparse(m.X.m)'*sparse(m.X.m)))
  #Hlinv = inv([ m.s.λtZt *m.s.λtZt'+eye(q) k12 ; k12' sparse(m.X.m)'*sparse(m.X.m)])
  spλtZt=sparse(m.s.λtZt)
  RZXa = spλtZt * m.X.m
  for j in 1:size(RZXa,2)
            permute!(view(RZXa,:,j),m.s.perm)
  end
  RZXa = Base.LinAlg.CHOLMOD.solve(m.s.L, RZXa, CHOLMOD_L)
  Hlinv2=inv(tril(hvcat((2,2),sparse(m.s.L),RZXa,RZXa',m.s.RX.UL')))       
  Hlinv= Hlinv2' * Hlinv2 
  rho = trace( vcat(spλtZt,sparse(m.X.m)')*vcat(spλtZt,sparse(m.X.m)')' *Hlinv)
  #REML
  #rho=(n-p-1.0)*(rho+1.0)/(n-p-2.0) + (p+1.0)/(n-p-2.0)
  #MLE
  rho=n/(n-p-2.0) * ((rho+1) -(rho-p)/(n-p))
  (-2.0*condll(m)+2.0*rho, rho) 
end

function objectivecon(m::LinearMixedModel)
    n,p,q,k = size(m); fn = float64(n - (m.REML ? p : 0))
    sig2 = std(m)[end][1]^2
    n*log(2*pi*sig2)+sumfdiff(Abs2Fun(), m.y, m.mu)/sig2
end

function TRCderHat(f::LinearMixedModel, fo::Formula, fr::AbstractDataFrame)
    rhohat=0.0
    n,p,q,k = size(f)
    
    for i in 1:n
      x=f.y[i]
      h = x == 0 ? sqrt(eps(Float64)) : sqrt(eps(Float64)) * x
      xph = x + h
      dx = xph - x
 
      newdf = deepcopy(fr)
      newdf[fo.lhs]=convert(DataArray{Float64,1},newdf[fo.lhs])
      newdf[fo.lhs][i]=xph
      f1= fit(lmm(fo,newdf,true)).mu[i]
      f0 = f.mu[i]
      rhohat += (f1 - f0) / dx
    end

    rhohat
end

function HAIC(m::LMMGeneral)
  n,p,q,k = size(m)
  k12= m.LambdatZt * sparse(m.X.m)  
  #Hlinv = inv(hvcat((2,2),m.LambdatZt *m.LambdatZt' +eye(q) ,k12,k12',sparse(m.X.m)'*sparse(m.X.m)))
  Hlinv = inv([ m.LambdatZt *m.LambdatZt'+eye(q) k12 ; k12' sparse(m.X.m)'*sparse(m.X.m)])
  vcat(m.LambdatZt,sparse(m.X.m)')' *Hlinv*vcat(m.LambdatZt,sparse(m.X.m)')
end

function TRCderHat2(f::LinearMixedModel, G::Array{Float64},
                    fo::Formula, fr::AbstractDataFrame)
    
    n,p,q,k = size(f)
    Hout = eye(n)

    for i in 1:n
      for j in 1:n
      x=f.y[j]
      # pick a small value for h
      h = x == 0 ? sqrt(eps(Float64)) : sqrt(eps(Float64)) * x
      # floating point arithmetic gymnastics
      xph = x + h
      dx = xph - x
 
      # evaluate f at x + h
      newdf = deepcopy(fr)
      #newdf[fo.lhs]=convert(DataArray{Float64,1},newdf[fo.lhs])
      newdf[fo.lhs][j]=xph
      f1= HAIC(fit(lmm(fo,newdf,true)))[i,j]
      # evaluate f at x
      f0 = G[i,j]
      # divide the difference by h
      Hout[i,j]= (f1 - f0) / dx
      #Hout[j,i]=Hout[i,j]
    end
    end
    ones(n)'*Hout*f.y
    #Hout
end

#????? ==> newdf = deepcopy(fr)??? 
function derHat(f::LinearMixedModel, idx::Int64, fo::Formula, fr::AbstractDataFrame)
    
    return function(x)
      
      # pick a small value for h
      h = x == 0 ? sqrt(eps(Float64)) : sqrt(eps(Float64)) * x
      # floating point arithmetic gymnastics
      xph = x + h
      dx = xph - x
 
      # evaluate f at x + h
      newdf = fr
      #newdf[fo.lhs]=convert(DataArray{Float64,1},newdf[fo.lhs])
      newdf[fo.lhs][idx]=xph
      f1= fit(lmm(fo,newdf,true)).mu[idx]
      # evaluate f at x
      f0 = f.mu[idx]
      # divide the difference by h
      (f1 - f0) / dx
    end

end