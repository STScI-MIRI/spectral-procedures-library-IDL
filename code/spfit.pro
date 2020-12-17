FUNCTION spfit,spa,spb,l0,l1,_extra=e

; 31 Dec 06 created
;
; spfit fits one spectrum (spa) to another (spb) over the wavelength
; range l0-l1 using calls to spavg and sptimes
;
; INPUT
;   spa - spectrum to fit
;   spb - spectrum to be fitted to
;   l0,l1 - wavelength range for multiplicative fit
; OUTPUT
;   returns spa after being shifted to fit spb

RETURN,sptimes(spa,spavg(spb,l0,l1,_extra=e)/spavg(spa,l0,l1,_extra=e))
END
