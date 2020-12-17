FUNCTION spdivide,a,b

; 31 Mar 04 now can take second argument as vector or scalar
; 28 Feb 04 removing bugs (GS)
; 31 Jan 04 brought into alignment with other spfits code (GS)
; 13 Nov 03 created by M.J. Russell
;
; Divides two spectra
; INPUT - a - numerator, as a spectral data array
;       - b - denominator, as a spectral data array or scalar
; OUTPUT - returns array a, except fluxes are ratios and errors have
;           been propagated, see Bevington (1969, p 62)
;
;  note that the only check to see if the arrays are compatible is to
;    compare the number of rows of data
;  ASSUMED:  both arrays have a flux and error column (as cols 1 and 2)
;            both arrays, if segmented, have the same structure

if (n_elements(b) eq 1) then begin ; scalar division
  outarray=a
  outarray[1:2,*] = a[1:2,*]/b
endif else begin                   ; vector division

  outarray=a
  len = n_elements(a[0,*])

  if (n_elements(b[0,*]) ne len) then begin
    print,$
      'Warning in spdivide.  Arrays have different sizes, returning numerator.'
  endif else begin

    num = a[1,*] & den = b[1,*]
    n_e = a[2,*] & d_e = b[2,*]

; replace zeroes with a small number to avoid dividing by zero

    index = where(den eq 0)
    if (max(index) gt -1) then den[index] = .00000001

    outarray[1,*] = num/den
    outarray[2,*] = sqrt(  outarray[1,*]^2*((n_e/num)^2 + (d_e/den)^2)  )

  endelse
endelse

return,outarray

END
