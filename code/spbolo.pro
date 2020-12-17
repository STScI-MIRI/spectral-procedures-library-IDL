FUNCTION spbolo,spin,la0,la1,lb0,lb1,T

; 20 Sep 13 don't add Wien tail if spectrum goes below 0.1 um
; 30 Apr 10 activated code to add in RJ tail - it was never added before!
;           also forced normalization of Wien tail to be zero or greater
;  4 Oct 05 created
;
; integrates a spectrum from zero to infinity
; extrapolates from zero to min wavelength with Planck f'n of temperature T
;   (T defaults to 3000 K)
;   actually pads spectrum with data starting at 1 A with 50 data per decade
; extrapolates from max wavelength to infinity with Rayleigh-Jeans tail
;   does this analytically:  Integrated flux = A / ( 3 l_max^3 )
;   where A normalizes RJ tail to spectrum
;
; INPUT
;   spin  - spectral data array (assume wavelength in um, flux in Jy)
;   la0,1 - wavelength range to fit Wien side of spectrum
;   lb0,1 - wavelength range to fit Rayleigh-Jeans tail
;   T     - temperature for Wien side (defaults to 3000 K)
; OUTPUT  - integrated flux in W m^-2

lcol=0 & fcol=1

; sort spin by wavelength

idx=sort(spin[lcol,*])
sp=spin[*,idx]


minlam=spin[0,0] ; true enough after sorting

if (minlam gt 0.10) then begin

; build wavelength grid for Wien-side spectrum
; grid runs from 1 Angstrom (1e-4 um) to 1 mm (1000 um)
;   with 10 steps per decade

  ncols=n_elements(sp[*,0])
  nrows=351
  new = fltarr(ncols,nrows)
  new[lcol,*] = 10 ^ (-4 + findgen(nrows) / 50)

; if T not set, set = 3000 K
; load Wien-side spectrum with a Planck f'n of temperature T (in F_nu units)

  if (n_elements(T) eq 0) then T=3000.0
  new=spplanck(new,T)

; normalize Wien-side spectrum to input spectrum

  laspot=mean([la0,la1])
  normal=spavg(sp,la0,la1)/spspot(new,laspot)
  if (normal lt 0) then normal=0 ; negative fluxes are not real!
  new=sptimes(new,normal)

; then cut Wien-side spectrum and concatenate input spectrum to it

  l_min=min(sp[lcol,*])
  new=spcut(new,0,l_min)
  new=[[new],[sp]]

endif else new=sp               ; do not include Wien side

; cut spectrum past lb1 and convert to F_lam units (W m^-2 um^-1)

new=spcut(new,0,lb1)
new[fcol,*]=new[fcol,*] * 3e-12 / new[lcol,*]^2
nrows=n_elements(new[lcol,*])

; integrate from lam = 0 to lb1

dlam=shift(reform(new[lcol,*]),-1)-reform(new[lcol,*])
dlam[nrows-1]=dlam[nrows-2]
flux_bol=total(reform(new[fcol,*])*dlam)

; find normalization for RJ tail by generating a fake spectrum
; use this to determine flux_rj, addition to flux_bol for RJ tail

rj=new
rj[fcol,*]=1.0/rj[lcol,*]^4
normal=spavg(new,lb0,lb1)/spavg(rj,lb0,lb1)
flux_rj = normal / (3 * lb1^3) ; analytical integral from l_max to infinity
flux_bol = flux_bol + flux_rj

RETURN,flux_bol
END
