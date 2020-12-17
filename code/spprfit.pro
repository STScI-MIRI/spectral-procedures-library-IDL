FUNCTION spprfit,sp,l0,l1

;  7 Jul 06 created 
;
; Given a spectrum, fit the Jura model Poynting Robertson 
; dust distribution
;
; INPUT
;   sp    - input spectral data array
;   l0,l1 - range of wavelengths to fit    
; OUTPUT  - returns fitted model

lcol=0 & fcol=1 & ecol=2 

; Compute model

spnew=sp
spnew[fcol,*]=spnew[lcol,*]
spnew[ecol,*]=0

; Fit model to spectrum

spnew = sptimes(spnew,spavg(sp,l0,l1)/spavg(spnew,l0,l1))

RETURN,spnew
END
