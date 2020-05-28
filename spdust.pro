FUNCTION spdust,sp,t,LAM=lam,ENGELKE=engelke

; 16 Dec 04 created
;
;  takes a spectral FITS array, fits and subtracts a continuum, 
;  returns the difference
; 
;  INPUT
;    sp      - spectral FITS array
;    t       - temperature of Planck or Engelke function
;    lam     - optional keyword to set wavelengths for continuum fitting
;              default is 7.60-8.60 (might should be 7.67-8.74...)
;    engelke - optional keyword to fit an Engelke f'n and not a Planck f'n

; check keywords, set defaults, generate continuum

if (keyword_set(lam) eq 0) then lam=[7.60,8.60]
if (keyword_set(engelke) eq 0) then csp=spplanck(sp,t) $
  else csp=spengelke(sp,t)

; normalize continuum spectrum to science spectrum and subtract it

factor = spavg(sp,lam[0],lam[1])/spavg(csp,lam[0],lam[1])
fitsp  = sptimes(csp,factor)
dust   = spadd(sp,fitsp,/minus)

RETURN,dust 

END
