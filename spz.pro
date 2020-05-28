FUNCTION spz,sp,z

; 16 Mar 05 repaired
;  8 Dec 04 created
;
;  given a spectral data array and a redshift, 
;  move the wavelength grid to the rest frame
;
;  INPUT
;    sp - input spectral data array
;    z  - redshift
;  OUTPUT - returns spectrum with wavelength shifted to the rest frame

; z = v/c = Dlam / lam = (lam-lam0)/lam0
; lam0 = lam/(1+z)

spnew=sp
spnew[0,*] = sp[0,*]/(1+z)

RETURN,spnew
END
