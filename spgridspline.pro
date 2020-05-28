FUNCTION spgridspline,sparray,lamarray

;  2 Nov 12 modified to handle sparray with only two input columns
; 21 Apr 07 guaranteeing at least 4 columns in output array
; 17 May 06 tidied up comments
; 27 Apr 06 major rewrite, no longer regridding order-by-order
; 30 Jan 04 created
;
; regrids a spectral data array using the IDL routine cspline
;   previous version performed regridding order-by-order
;     this required a match between lamarray and sparray
;   that is no longer necessary - updated version sorts sparray by wavelength
;     must be passed to cspline.pro with monotonically changing grid
;       now passed to cspline.pro as one large array
;     lamarray no longer has to be monotonically changing
;
; INPUT  - sparray  - spectral data array to be regridded
;          lamarray - spectral data array with new wavelength grid
;                     if this array has a col 3, this is assumed to be order
;                     otherwise, col 1 is assumed to be the order column
; OUTPUT - returns a new spectra data array with wavelength,
;          regridded fluxes and errors, output order

; set up, check dimensions, create dummy vector (in case it's needed)

lcol=0 & fcol=1 & ecol=2 & ocol=3

ncol    = n_elements(sparray[*,0])
nrow    = n_elements(sparray[0,*])
nlamcol = n_elements(lamarray[*,0])
nlamrow = n_elements(lamarray[0,*])

dummy   = fltarr(nrow)
if (nlamcol ge 4) then olcol=3 else olcol=1

; copy sparray to gridarray and sort by wavelength

testarray = sparray
idx = sort(sparray[lcol,*])
testarray = sparray[*,idx]
  
old_lam  = reform(testarray[lcol,*])
old_flux = reform(testarray[fcol,*])
if (ncol gt ecol) then old_err =reform(testarray[ecol,*]) else old_err=dummy
new_lam  = reform(lamarray [lcol,*])

new_flux = cspline(old_lam,old_flux,new_lam)
new_err  = cspline(old_lam,old_err ,new_lam)

outarr = fltarr(4,nlamrow)
outarr[lcol,*] = lamarray[lcol,*]
outarr[fcol,*] = new_flux
outarr[ecol,*] = new_err
outarr[ocol,*] = lamarray[olcol,*]

return,outarr
END
