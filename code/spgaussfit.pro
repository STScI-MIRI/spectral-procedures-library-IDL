FUNCTION spgaussfit,spec,lam0,lam1,NTERMS=nterms,_extra=e,PLOT=plot

; 13 Jul 12 modified wavelengths, so that we never change input lam0 and lam1
;  3 Jan 06 added nterms parameter, check that we have more data than terms
; 19 Apr 05 cleaning up
; 15 Apr 05 added over flag
; 14 Apr 05 created
;
; fits a Gaussian to a spectral data array from lam0 to lam1
; returns the fitting parameters [l0,sigma,amplitude]
; _e allows nterms to be passed to gaussfit
;
; INPUT
;   spec       - spectral data array
;   lam0,lam1  - wavelength range to consider (optional, but recommended)
;   plot       - optional keyword to plot results
;   nterms     - optional keyword for gaussfit, defaults to 3
;   _extra - used to pass keywords on to gaussfit, spplot, and oplot
;          - use /over to overplot on a plot already generated before spgauss
; OUTPUT
;   returns an array with the following:
;     if nterms=3, amplitude,lam_0,sigma
;     if nterms=4, amplitude,lam_0,sigma,flux_offset
;     if nterms=5, amplitude,lam_0,sigma,y_int,slope

lcol=0 & fcol=1

; set nterms keyword if necessary

if (n_elements(nterms) eq 0) then nterms=3

; check wavelengths and set if necessary

if (n_elements(lam1) eq 0) then begin
  l0=min(spec[lcol,*])
  l1=max(spec[lcol,*])
endif else begin
  l0=lam0 & l1=lam1
endelse
idx=where(spec[lcol,*] ge l0 and spec[lcol,*] le l1)

; fit the gaussian over the indexed wavelength range

if (idx[0] gt -1 and n_elements(idx) gt nterms) then begin

  lam=double(reform(spec[lcol,idx]))
  flux=double(reform(spec[fcol,idx]))
  fit=gaussfit(lam,flux,coeff,nterms=nterms,_extra=e)
  if (n_elements(plot) gt 0) then begin
    if (n_elements(over) eq 0) then spplot,spec,xran=[l0,l1],_extra=e $
      else spplot,spec,_extra=e,psym=-1
    oplot,lam,fit,_extra=e,th=2
  endif
  retarr=coeff

; or return an error message

endif else begin
  print,'Error in spgaussfit'
  if (idx[0] lt 0) then print,'No data in range ',l0,' - ',l1
  if (n_elements(idx) lt nterms) then print,'Not enough data to fit'
  retarr=fltarr(nterms)
endelse

return,retarr
END
