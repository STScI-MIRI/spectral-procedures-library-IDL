FUNCTION spfspot,spfile,lam,_extra=e

; 29 Oct 04 created
;
; reads in spfile, uses spspot at wavelength lam to find flux
; prints flux and ends

sp=readfits(spfile,/silent)
flux=spspot(sp,lam,_extra=e)

return,flux

END
