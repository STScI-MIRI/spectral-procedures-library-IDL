FUNCTION spstar,inspec,T,lam0,lam1,PLANCK=planck

; 17 Sep 10 created
;
; spstar generates a stellar continuum of temperature T fitted to inspec
;   in the interval lam0-lam1
; keyword PLANCK changes the continuum from a default Engelke function

if (keyword_set(planck) eq 0) then eflag=1 else eflag=0

if (eflag eq 1) then cont=spengelke(inspec,T) else cont=spplanck(inspec,T)

cont=sptimes(cont,spavg(inspec,lam0,lam1)/spavg(cont,lam0,lam1))

RETURN,cont
END
