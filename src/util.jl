 function checkzerovar(m::LinearMixedModel)
    ff = deepcopy(m.f)
    # set fixed term=true, reterm=false in rhs of formula

    fixterms=map(x-> !(Meta.isexpr(x,:call) && x.args[1] == :|), m.mf.terms.terms)
    
    retrms = filter(x->Meta.isexpr(x,:call) && x.args[1] == :|, m.mf.terms.terms)
    eachtrms=[aa.args[2] for aa in retrms]
    grps = Any[t.args[3] for t in retrms]
    ugrps = unique(grps)

    #check the numer of terms in each group
    trmsleng=zeros(length(eachtrms))

    for k in 1:length(eachtrms)
      if isa(eachtrms[k],Int64)  
        trmsleng[k]=1
      elseif  isa(eachtrms[k],Symbol)
          trmsleng[k]=2 
      elseif  eachtrms[k].args[1] == :+ && eachtrms[k].args[2] == 0
        trmsleng[k]=length(eachtrms[k].args)-2
      else
        trmsleng[k]=length(eachtrms[k].args)
      end 
    end


    #find λ and their lowe bounds , compare them, 
    # set true if boundary conditions are met or set false otherwise                
    est=map(MixedModels.θ,m.λ)
    lowbd=map(MixedModels.lower,m.λ)
    

    if all(trmsleng .== 1.0)

      if length(ugrps) == length(grps) 
        remainterms1=map(all,Any[l .> k  for (l,k) in zip(est,lowbd)])
    
      else # length(ugrps) < length(grps) 
        aa=[l .> k  for (l,k) in zip(est,lowbd)]
        remainterms1=aa[1]
        for i in 2:length(aa)
         remainterms1=[remainterms1,aa[i]]
        end      
      end

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
      return(ff)

    else

      ffmt=deepcopy(m.mf.terms.terms)
      numfixterm= sum(fixterms .== true)

      if length(ugrps) == length(grps)  
         
          aa=[l .> k  for (l,k) in zip(est,lowbd)]
         
          for i in 1:length(aa)
              if length(aa[i])==1 
                
                if aa[i][1]==true
                  fixterms[numfixterm+i]=true
                else #  aa[i][1]==false
                  fixterms[numfixterm+i]=false
                end
   
         #  especially for model m5 
              elseif length(aa[i])==3 
                   
                   if aa[i][1] == true && aa[i][3] == true 
                     fixterms[numfixterm+i]=true
                   elseif aa[i][1] == false && aa[i][3] == true 
                     fixterms[numfixterm+i]=true
                     ffmt[2]= :(0+time | id)
                   elseif aa[i][1] == true && aa[i][3] == false 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :(1 | id)
                   else 
                     fixterms[numfixterm+i]=false 
                   end 
         #  especially for model m7 
   
              elseif length(aa[i])==6 
                   
                   if aa[i][1] == true && aa[i][4] == true  && aa[i][6] == true 
                     fixterms[numfixterm+i]=true
                   
                   elseif aa[i][1] == false && aa[i][4] == true  && aa[i][6] == true 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :((0+time + time2) | id)
                   elseif aa[i][1] == true && aa[i][4] == false  && aa[i][6] == true 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :((1+ time2) | id)
                   elseif aa[i][1] == true && aa[i][4] == true  && aa[i][6] == false 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :((1+ time) | id)
                   
                   elseif aa[i][1] == false && aa[i][4] == false  && aa[i][6] == true 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :((0+ time2) | id)
                   elseif aa[i][1] == true && aa[i][4] == false  && aa[i][6] == false 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :(1 | id)
                   elseif aa[i][1] == false && aa[i][4] == true  && aa[i][6] == false 
                     fixterms[numfixterm+i]=true
                     ffmt[2] = :((0+ time) | id)
       
                    else 
                      fixterms[numfixterm+i]=false 
                    end 
   
               end
           end
        end

        ffm=ffmt[convert(Array{Bool,1},fixterms)]
      
        a=Expr(:call,tuple ( [:+,ffm]...)...)
        ff=eval(y~$a)
    

      return(ff)
    
    end

  
end
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


