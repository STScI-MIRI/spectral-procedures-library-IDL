FUNCTION spcontp,sp,l0,l1,temp

; 27 Apr 05 created
;
; spcontp fits a Planck function to a spectrum from l0 to l1
; if temp not supplied, temp=10000 K

if (keyword_set(temp) eq 0) then temp=10000.0

cont=spplanck(sp,temp)
cont=sptimes(cont,spavg(sp,l0,l1)/spavg(cont,l0,l1))

RETURN,cont
END
