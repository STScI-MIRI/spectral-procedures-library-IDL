FUNCTION spadd,a,b,minus=MINUS,noerror=NOERROR

;  4 Aug 06 added noerror keyword to turn off error propagation 
; 26 Apr 06 modified to work if a has error column but b does not
; 23 Mar 05 modified to work if there is no error column in a
; 31 Mar 04 created from spdivide.pro
;
; adds or subtracts two spectra or one spectrum and a scalar
; INPUT - a - spectral data array
;       - b - spectral data array or scalar
; OUTPUT - returns array a, except fluxes are sums or differences and 
;          errors have been propagated, see Bevington (1969, p 62)
;
;  note that the only check to see if the arrays are compatible is to
;    compare the number of rows of data
;  ASSUMED:  both arrays have a flux and error column (as cols 1 and 2)
;            both arrays, if segmented, have the same structure

lcol=0 & fcol=1 & ecol=2 & ocol=3

if (keyword_set(minus) eq 0) then addsign=1 else addsign=-1

; 
if (n_elements(b) eq 1) then begin ; scalar addition

  outarray=a
  outarray[fcol,*] = a[fcol,*]+addsign*b

endif else begin                   ; vector addition

  outarray=a
  len = n_elements(a[lcol,*])

  if (n_elements(b[lcol,*]) ne len) then begin
    print,'Warning in spadd.  Arrays have different sizes, returning first argument.'
  endif else begin

    v1 = a[fcol,*] 
    v2 = b[fcol,*]
    outarray[fcol,*] = v1 + addsign*v2

;   propagate errors if error col present in both arrays 
;     and noerror keyword not set

   if (n_elements(a[*,0]) gt 2 and n_elements(b[*,0]) gt 2 $
     and keyword_set(noerror) eq 0) then begin
     e1 = a[ecol,*] 
     e2 = b[ecol,*]
     outarray[ecol,*] = sqrt(e1*e1 + e2*e2)
    endif

  endelse
endelse

return,outarray

END
