type IC <: InfoCrit
   likelihood ::Float64
   df         ::Float64
   method     ::ASCIIString
   crterm     ::Vector{Float64}
   value      ::Float64 
end

