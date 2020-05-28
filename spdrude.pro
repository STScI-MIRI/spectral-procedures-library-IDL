FUNCTION spdrude,lamarray,amp,lam0,fwhm

; 31 Jul 12 created
;
; generates a Drude profile 
; reference:  Draine & Malhotra (1993, ApJ, 414, 632)
; INPUT
;   lamarray  - a spectral data array with wavelength in col 0 (in um)
;   amp       - amplitude in arbitrary units
;   lam0      - central wavelength in um
;   fwhm      - FWHM (aka gamma) in arbitrary units

; check lamarray to determine which column is order

ncol=n_elements(lamarray[*,0])
nrow=n_elements(lamarray[0,*])
lcol=0
if (ncol lt 4) then ocol=1 else ocol=3

c=1 ; should be the speed of light (2.9979e14 um/second)
    ; but set=1 to simplify normalization

outarray=fltarr(4,nrow)
lam=reform(lamarray[lcol,*])

drudedata = amp * c * fwhm^2 / ( fwhm^2 + (1/lam - lam/lam0^2)^2 )

outarray[0,*]=lam
outarray[1,*]=drudedata
outarray[3,*]=reform(lamarray[ocol,*])

RETURN,outarray
END
