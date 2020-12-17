FUNCTION spregrid,sparray,lamarray

; 29 Nov 12 now just copying old to new if the wavelength grids are the same
;  2 Nov 12 made all for loop counters long and sparray with only 2 columns
; 13 Jul 09 renamed spregrid
;  8 Jul 09 renamed spshift, modified sm_ext_opt_spgridspline by G. Sloan
; 10 Jun 07 modified by D. Whelan for use within SMART
; 30 Jan 04 created as spgridspline by G. Sloan
;
; regrids a spectral data array using splines based on the integral of the flux
;   spgridspline used cspline; this routine uses calls to interpol
;   previous version performed regridding order-by-order
;     this required a match between lamarray and sparray
; lamarray does not have to be monotonic
;
; INPUT  - sparray  - spectral data array to be regridded
;          lamarray - spectral data array with new wavelength grid
;                     if this array has a col 3, this is assumed to be order
;                     otherwise, col 1 is assumed to be the order column
; OUTPUT - returns a new spectra data array with wavelength,
;          regridded fluxes and errors, output order

; define columns, determine dimensions, create dummy array (if needed)

lcol=0 & fcol=1 & ecol=2 & ocol=3

ncol    = n_elements(sparray[*,0])
nrow    = n_elements(sparray[0,*])
len     = n_elements(lamarray[lcol,*])
if (n_elements(lamarray[*,0]) ge 4) then olcol = 3 else olcol = 1
dummy   = fltarr(nrow)
new_lam = reform(lamarray[lcol,*])

; copy sparray into sortarray, sort by wavelength, and load old vectors

sortarray  = sparray
sidx       = sort(sparray[lcol,*])
sortarray  = sparray[*,sidx]
old_lam    = reform(sortarray[lcol,*])
old_flux   = reform(sortarray[fcol,*])
if (ncol gt ecol) then old_err=reform(sortarray[ecol,*]) else old_err=dummy

; find the average wavelength spacing

avg_bin=0.0
for i=long(1),nrow-1 do avg_bin = avg_bin + (old_lam[i]-old_lam[i-1])
avg_bin = avg_bin / (nrow-1)

; drop NaN and Inf data from old vectors

idx = where(finite(old_lam) and finite(old_flux))
old_lam  = old_lam[idx]
old_flux = old_flux[idx]
old_err  = old_err[idx]
old_len  = n_elements(idx)

; check wavelength grids - regrid only if they differ
;   otherwise just copy old to new

diff=abs(old_lam-new_lam)
if (min(diff) gt 0 or max(diff) gt 0) then begin

; generate an integral spectrum using tsum
;   trapezoidal summation, in /usr/local/rsi/idl_library/astrolib/pro/math/
; can't use int_tabulated - not adequate for oscillatory-non smooth data
; interpolate the integral spectrum to the new wavelength grid
; and take the derivative

  int_old_flux = dblarr(n_elements(old_flux))
  for i=long(1),old_len-1 do int_old_flux[i] = tsum(old_lam[0:i],old_flux[0:i])
  int_new_flux = interpol(int_old_flux,old_lam,new_lam,/lsquadratic)
  new_flux     = deriv(new_lam,int_new_flux)

; repeat for the errors

  int_old_err = fltarr(old_len)
  for i=long(1),old_len-1 do int_old_err[i] = tsum(old_lam[0:i],old_err[0:i])
  int_new_err = interpol(int_old_err,old_lam,new_lam,/lsquadratic)
  new_err     = deriv(new_lam,int_new_err)

endif else begin

  new_flux = old_flux
  new_err  = old_err

endelse

; load output array, 5 columns (last used for flags and currently blank)

outarr = fltarr(5,len)
outarr[lcol,*] = new_lam
outarr[fcol,*] = new_flux
outarr[ecol,*] = new_err
outarr[ocol,*] = lamarray[olcol,*]

RETURN,outarr
END
