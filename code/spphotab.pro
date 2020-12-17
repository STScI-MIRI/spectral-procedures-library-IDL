FUNCTION spphotab,spin,filter,zp,FNU=fnu

; 12 Feb 14 created by modifying spphot.pro
;           now following AB mag formulation of Bessell (2012) for optical
; 25 Nov 12 changed loop counters to long to handle big spectral FITS files
;  2 Aug 07 implemented zeromag keyword - flux at zero magnitude
;  6 Feb 06 minor changes to comments
; 17 Oct 05 created as spphot.pro

; spphotopt multiplies a spectral data array by a filter response function
;   and integrates the product sp * filter / nu or sp * filter / lam

; NOTE - no error propagation 

; INPUT
; spin    - spectral data array - expects lambda (um) and F_nu (Jy)
; filter  - filter response function, col 0 = wavelength, col 1 = response
; zp      - zero point of filter in question (cf Bessell 2012)
; fnu     - keyword to require integration in frequency space (vs. wavelength)
; NOTE    - fnu NOT implemented

; algorithm:
; 0) - initialization
; 1) setup
;    - check keywords and initialize
;    - sort spectral data array by frequency 
;    - convert both to F_lam units
;    - regrid filter function to spectral data array
;    - compute multiply the two spectra by the filter function
;    - generate delta lambda vector (dellam)
; 2) integrate F_lam dlam and filter function - units of W m^-2
; 3) calibrate (take -2.5*log, subtact -48.60, add zero-point mag.)

; (0) check keywords and initialize

if (keyword_set(zp) eq 0) then zp=0.0

lcol=0 & fcol=1 & ecol=2 & ocol=3

; (1) setup input spectrum and filter function

sp  = spsort(spin,/nodupe) ; sort input spectrum
lam = sp[lcol,*]
len = n_elements(lam)

; convert input spectrum to F_lam units - W m^-2 um^-1

sp[fcol,*]   = sp[fcol,*]   * 3e-12 / (lam^2)
sp[ecol,*]   = sp[ecol,*]   * 3e-12 / (lam^2)

; regrid filter f'n and set to zero outside given range

spf = spgridspline(filter,sp)
idx=where( lam lt min(filter[lcol,*]) or lam gt max(filter[lcol,*]) )

if (max(idx) gt -1) then spf[fcol,idx] = 0.0

spx = sptimes(sp,spf)      ; multiply input spectrum by filter f'n 
dellam=abs(reform(shift(lam,1)-lam)) ; delta lambda vector
dellam[0]=dellam[1]                                ; must replace first pixel

; (2) integrate spectrum, filter function, Rayleigh-Jeans tail

sumsp=0.0 & sumff=0.0
for i=long(0),len-1 do begin
  sumsp = sumsp + spx[fcol,i] * dellam[i] * lam[i]
  sumff = sumff + spf[fcol,i] * dellam[i] / lam[i]
endfor

; (3) calibrate the integral and divide by sumff
;     normalize - divide by c, multiply by 0.1 (to go to erg cm^-2 s^-1 A^-1)

sumsp = sumsp * 0.1 / (3e10 * sumff)

out = -2.5*alog10(sumsp) - 48.60 - zp

RETURN,out
END
