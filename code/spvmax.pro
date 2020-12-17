FUNCTION spvmax,sp1,sp2
;
; 22 Dec 05 modified for long spectral arrays
; 16 Dec 05 created 
;
; spvmax returns the maximum value of two spectra at each wavelength.
; for the error, it returns the error corresponding to the max value
;
; INPUT
;   sp1 - input spectral array 1
;   sp2 - input spectral array 2
; OUTPUT - the new maximum spectrum (as a spectral array)

lcol=0 & fcol=1 & ecol=2

spnew=sp1

if (n_elements(sp1) eq n_elements(sp2)) then begin
  for i=long(0),n_elements(sp1[lcol,*])-1 do $
    if (sp2[fcol,i] gt sp1[fcol,i]) then spnew[*,i]=sp2[*,i]
endif else begin
  print,'Error in spvmax, spectra have different sizes, returning the first' 
endelse

RETURN,spnew
END
