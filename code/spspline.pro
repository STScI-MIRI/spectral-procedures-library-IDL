FUNCTION spspline,spin,lambda,spotarray,SMOOTHBOX=smoothbox,DIAG=diag,SUBTRACT=subtract,_extra=e

; 29 Dec 09 added /subtact keyword
; 31 Jul 07 added spotarray variable to return spline wavelengths and fluxes
; 29 Apr 06 created
;
; given an input spectrum and a series of wavelengths,
;   find the flux at those wavelengths using spspot
;   and fit a spline through the fluxes to generate a smooth output spectrum
;
; INPUT
;   spin      - the input spectral data array
;   lambda    - the series of wavelengths where the spline will be fit
;               these wavelengths cannot be outside the range of spin
;   spotarray - to return wavelengths and fitted fluxes (as spectral data array)
;   smoothbox - keyword to turn on smoothing
;   subtract  - keyword to return the spectrum minus the spline
;   diag      - keyword to turn on plotting

lcol=0 & fcol=1

; load spsmoo with contents of spin, smoothing of requested

if (keyword_set(smoothbox) gt 0) then spsmoo=spsmooth(spin,smoothbox,/noorder) $
  else spsmoo=spin

; first set up lam and flux vectors to pass to cspline

nlam=n_elements(lambda)
count=0
for i=0,nlam-1 do begin
  lam0=lambda[i]
  if (lam0 ge min(spsmoo[lcol,*]) and lam0 le max(spsmoo[lcol,*])) then begin
    flux0=spspot(spsmoo,lam0)
    if (count eq 0) then begin
      lam=lam0       & flux=flux0
    endif else begin
      lam=[lam,lam0] & flux=[flux,flux0]
    endelse
    count=count+1
  endif
endfor
nlam=count

; create spotarray - spectral data array with wavelengths and fitted fluxes

spotarray=fltarr(4,n_elements(lam))
spotarray[0,*]=lam
spotarray[1,*]=flux

; next set the output grid to be the wavelength grid from spsmoo
; then call cspline to generate a new set of fluxes
; finally load these fluxes into spout and return

if (keyword_set(diag) eq 1) then begin
  spplot,spin
  spplot,spsmoo,/over
endif

lamgrid=reform(spsmoo[lcol,*])
newflux=cspline(lam,flux,lamgrid)
spout=spsmoo
spout[fcol,*]=newflux

if (keyword_set(diag) eq 1) then begin
  spplot,spotarray,psym=4,symsize=2,th=2
  spplot,spout,/over,th=2
endif

if (keyword_set(subtract) ne 0) then spout=spadd(spsmoo,spout,/minus)

RETURN,spout
END
