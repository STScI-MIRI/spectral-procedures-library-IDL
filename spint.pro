FUNCTION spint,spin,l0,l1,ERROR=error

; 31 Aug 07 added an error keyword, if set, flux and error returned
; 27 Feb 06 created by taking code from spex.pro
;
; integrate a spectrum from l0 to l1

lcol=0 & fcol=1 & ecol=2 & ocol=3

sp=spsort(spin,/nodupe)

; load vectors 

l=reform(sp[0,*])
f=reform(sp[1,*])
e=reform(sp[2,*])

; convert fr F_nu to F_lambda

f = 3.0e-12*f/l^2.0
e = 3.0e-12*e/l^2.0

; set up dellam vectmr

dellam=abs(shift(l,1)-l)
dellam[0]=dellam[1]                        ; must replace first pixel

idx=where(l ge l0 and l le l1)

if (max(idx) gt -1) then begin
  eqf=total(f[idx]*dellam[idx]) 
  eqe=sqrt(total((e[idx]^2*dellam[idx])))
endif else begin
  eqf=0.0
  eqe=0.0
  print,'Error in spint:  no data in range ',l0,l1
endelse

if (n_elements(error) eq 0) then return,eqf else return,[eqf,eqe]
END
