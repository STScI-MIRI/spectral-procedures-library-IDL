FUNCTION sptimes,a,b,ORDER=order

; 12 Apr 05 modified to run in the event that only one array has errors
; 29 Dec 04 added keyword ORDER to limit correction to a single order
; 27 May 04 modified to handle arrays with no error column
; 31 Mar 04 created from spdivide.pro
;
; multiplies two spectra or one spectrum by a scalar
; INPUT - a - spectral data array
;       - b - spectral data array or scalar
; OUTPUT - returns array a, except fluxes are products and errors have
;           been propagated, see Bevington (1969, p 62)
;
;  note that the only check to see if the arrays are compatible is to
;    compare the number of rows of data
;  ASSUMED:  both arrays have a flux and error column (as cols 1 and 2)
;            both arrays, if segmented, have the same structure

lcol=0 & fcol=1 & ecol=2 & ocol=3

sz=size(a)
outarray=a

if (n_elements(order) ne 0) then idx=where(a[ocol,*] eq order) $
  else idx=findgen(sz[2])

if (max(idx) gt -1) then begin

  if (n_elements(b) eq 1) then begin                ; scalar division
    outarray[fcol,idx] = a[fcol,idx]*b
    if (sz[1] gt 2) then outarray[ecol,idx] = a[ecol,idx]*b
  endif else begin                                  ; vector division

    len = n_elements(a[lcol,idx])
    if (n_elements(b[lcol,idx]) ne len) then begin
      print,'Warning in sptimes.  Arrays have different sizes; no changes.'
    endif else begin

      v1=a[fcol,idx] & v2=b[fcol,idx]
      outarray[fcol,*] = v1*v2
      if (sz[1] gt 2) then begin
        if (n_elements(b[*,0]) gt 2) then begin
          e1=a[ecol,idx] & e2=b[ecol,idx] 
          outarray[ecol,*] = sqrt( e1*e1*v2*v2 + e2*e2*v1*v1 )
	endif else outarray[ecol,*]=a[ecol,idx]
      endif
    endelse
  endelse

endif else begin
  print,'Warning in sptimes.  Order ',order,' not in array; no changes.'
endelse

RETURN,outarray
END
