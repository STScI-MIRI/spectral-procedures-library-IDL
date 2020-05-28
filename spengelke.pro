FUNCTION spengelke,sparr,teff,FIT=fit,NONORM=nonorm,_extra=e
;
;  created long, long ago as planckn.pro
;
;   6 Jan 11 added _extra=e to allow passing of keywords to spfit
;   2 Jan 11 updated with nonorm keyword as in spplanck.pro
;   7 May 08 added fitting keyword to fit Planck f'n to sparr
;  12 Dec 07 if input spectral array 1-D, still returns 4 columns
;   4 Oct 07 zero error column in sp array, if it exists
;  13 Dec 04 modified from spplanck.pro (KEK)
;   6 Dec 04 modified from planckn.pro
;   6 Dec 04 changed loop counter to long to handle really large spectra
;  31 Dec 01 modified to allow lambda to be a column of a 2-D array
;
;  spengelke function returns the Engelke function in f_nu units
;    
;  INPUT
;     sparr  - spectral data array with col 0 = wavelength in um
;     teff   - effective blackbody temperature in K
;     fit    - optional fitting wavelengths for normalization of output to sparr
;     nonorm - to return physically meaningful data
;  OUTPUT - returns spectral data array with F_nu in col 0
;

lcol=0 & fcol=1 & ecol=2 & ocol=3

sz=size(sparr)
if (sz[0] ne 1) then spnew=sparr else begin
  spnew=fltarr(4,n_elements(sparr))
  spnew[0,*]=sparr
endelse
len=n_elements(spnew[lcol,*])
lam=reform(spnew[lcol,*])
flux=reform(spnew[fcol,*])
if (n_elements(spnew[*,0]) gt ecol) then spnew[ecol,*]=0.0

hc = 1.986e-16
c  = 2.9979e10
k  = 1.381e-16


for i=long(0),len-1 do begin
  l = lam[i]*1.0e-4 ; translate to cm
  lb= lam[i]        ; except Engelke (1992) defined constants with [l]=micron
  tb = 0.738 * teff * (1+ (79450./(lb*teff)))^0.182 ; modified temperature
  argument = double(hc/(l*k*tb))
  feng  = 2.0*hc/(l*l*l * (exp(argument)-1.0) ) ; f_nu units
  flux[i] = feng
endfor

if (keyword_set(nonorm) eq 0) then flux = flux / max(flux)
spnew[fcol,*]=reform(flux)

if (n_elements(fit) eq 2) then spnew=spfit(spnew,sparr,fit[0],fit[1],_extra=e)

RETURN,spnew
END
