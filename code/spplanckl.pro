FUNCTION spplanckl,sparr,T,NONORM=nonorm
;
;   7 Dec 10 added the nonorm keyword to keep everything in CGS units
;  28 Aug 06 created by modifying spplanck.pro
;
;  spplanckl function returns the Planck function in f_lam units
;
;  INPUT
;     sparr - spectral data array with col 0 = wavelength in um
;     T     - blackbody temperature in K
;  OUTPUT - returns spectral data array with F_nu in col 0
;

len=n_elements(sparr[0,*])
flux=sparr[1,*]
spnew=sparr

hc = 1.986e-16
c  = 2.9979e10
k  = 1.381e-16

for i=long(0),len-1 do begin
  l = sparr[0,i]*1.0e-4
  argument = double(hc/(l*k*T))
  B  = 2.0*hc*c/(l^5 * (exp(argument)-1.0) ) ; CGS units
  flux[i] = B
endfor

if (keyword_set(nonorm) eq 0) then flux = flux / max(flux)
spnew[1,*]=flux

RETURN,spnew
END
