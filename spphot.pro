FUNCTION spphot,spin,filter,lamf,ZEROMAG=zeromag,RJ=rj

; 11 Dec 15 returning 99.999
; 10 Sep 15 finally aligned normalization with IRAS and IRAC conventions
;           added /rj keyword to keep old normalization to a Rayleigh-Jeans tail
; 25 Nov 12 changed loop counters to long to handle big spectral FITS files
;  2 Aug 07 implemented zeromag keyword - flux at zero magnitude
;  6 Feb 06 minor changes to comments
; 17 Oct 05 created

; spphot multiplies a spectral data array by a filter response function
;   and integrates the product
; NOTE - no error propagation 

; INPUT
; spin    - spectral data array
; filter  - filter response function, col 0 = wavelength, col 1 = response
; lamf    - effective wavelength of the filter
; zeromag - flux at zero magnitude (Jy) - if set, returns magnitude, not flux
; rj      - keyword to set RJ normalization (vs. F_nu ~ lambda)

; algorithm:
; 1) setup
;    - sort spectral data array by wavelength 
;    - generate a Rayleigh-Jeans tail (set to 1.0 Jy at wavelength lamf)
;    - convert both to F_lam units
;    - regrid filter function to spectral data array
;    - multiply the two spectra by the filter function
;    - generate delta lambda vector (dellam)
; 2) setup calibration spectrum
;    - a Rayleigh-Jeans tail scaled to 1 Jy at lamf
; 3) integrate F_lam dlam and filter function - units of W m^-2
; 4) calibrate 
;    - convert observed integrated flux to Jy at lamf
; 5) convert to magnitude if zeromag set

  lcol=0 & fcol=1 & ecol=2 & ocol=3

; (1) and (2) setup input and calibration spectrum

  sp  = spsort(spin,/nodupe) ; sort input spectrum

; generate the normalization spectrum (called sprj for historical reasons)
; if rj set, Rayleigh-Jeans tail normalized to 1.0 Jy at wavelength lamf
; otherwise, F_nu ~ lambda

  sprj = sp
  if (keyword_set(rj) ne 0) then sprj[fcol,*]=1.0/(sprj[lcol,*]^2) $
    else sprj[fcol,*] = sprj[lcol,*]
  sprj=sptimes(sprj,1.0/spspot(sprj,lamf))

; convert input spectrum and normalization spectrum to F_lam units

  sp[fcol,*]   = sp[fcol,*]   * 3e-12 / (sp[lcol,*]^2)
  sprj[fcol,*] = sprj[fcol,*] * 3e-12 / (sp[lcol,*]^2)
  sp[ecol,*]   = sp[ecol,*]   * 3e-12 / (sp[lcol,*]^2)
  sprj[ecol,*] = 0.0

; regrid filter f'n and set to zero outside given range

  spf = sptem(sp,filter)     
  idx=where( spf[lcol,*] lt min(filter[lcol,*]) or $
             spf[lcol,*] gt max(filter[lcol,*]) )

  if (max(idx) gt -1) then spf[fcol,idx] = 0.0

  spx = sptimes(sp,spf)      ; multiply input spectrum by filter f'n 
  spy = sptimes(sprj,spf)    ; ditto for normalization spectrum
  dellam=abs(reform(shift(sp[lcol,*],1)-sp[lcol,*])) ; delta lambda vector
  dellam[0]=dellam[1]                                ; must replace first pixel

; (3) integrate spectrum, filter function, normalization spectrum

  len=n_elements(sp[lcol,*])
  sumsp=0.0 & sumff=0.0 & sumrj=0.0
  for i=long(0),len-1 do begin
    sumsp = sumsp + spx[fcol,i]*dellam[i]
    sumff = sumff + spf[fcol,i]*dellam[i]
    sumrj = sumrj + spy[fcol,i]*dellam[i]
  endfor

  sumsp=sumsp/sumff   ; normalize integral by dividing by integral of filter f'n
  sumrj=sumrj/sumff   ; ditto for integral of normalization spectrum
                      ; (these will cancel out in the next step...)

; (4) calibrate integral back to Jy by dividing by integral of norm. spectrum

  sumsp=sumsp/sumrj

; (5) convert to magnitude if zeromag set, set mag=99 if flux negative

  if (keyword_set(zeromag) ne 0) then begin
    if (sumsp gt 0) then sumsp = -2.5*alog10(sumsp/zeromag) else sumsp=99.999
  endif

RETURN,sumsp
END
