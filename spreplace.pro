FUNCTION spreplace,sparray,lambda,newflux,ORDER=order,ERROR=error

; 17 Aug 04 created
;
; spreplace replaces the flux for the wavelength element closest to
;   the input lambda
; if lambda is midway between two wavelength elements, then the latter one is
;   replaced
; if lambda outside spectral range, an endpoint will be replaced if it is
;   closer than one wavelength interval to lambda
; INPUT
;   sparray - spectral data array
;   lambda  - wavelength element to replace
;   newflux - new flux value for wavelength element
;   order   - optional, order to search for lambda
;   error   - optional, replaces error column instead of flux column
; NOTE spectral data in each order must be monotonic by wavelength

newarray=sparray

; set up column and order information

lcol=0
fcol=1
ecol=2
ocol=3
if (keyword_set(error) eq 0) then rcol=fcol else rcol=fcol
minorder=min(sparray[ocol,*])
maxorder=max(sparray[ocol,*])

if (keyword_set(order) eq 1) then begin
  minorder=order
  maxorder=order
endif

for m=minorder,maxorder do begin
  idx=where(sparray[ocol,*] eq m)
  if (max(idx) ge 0) then begin

    lthan=where(sparray[lcol,idx] lt lambda)
    gthan=where(sparray[lcol,idx] gt lambda)

    if (max(lthan) ge 0 and max(gthan) ge 0) then begin 

;     lambda inside order, replace flux or error value

      if (max(lthan) lt min(gthan)) then begin ; wavelength increasing
        i0=max(lthan) & i1=min(gthan)
      endif else begin                         ; wavelength decreasing
        i0=min(lthan) & i1=max(gthan)
      endelse

;     replace the closer (or latter) value

      if (abs(sparray[lcol,i0]-lambda) lt abs(sparray[lcol,i1]-lambda)) then $
         newarray[rcol,i0]=newflux else newarray[rcol,i1]=newflux

    endif else begin

;     lambda outside order
;     replace one endpoint if lambda closer to it than the adjacent wavelength

      i0=min(idx) & i1=max(idx)
      testlam=sparray[lcol,i0]
      if (abs(testlam-lambda) lt abs(testlam-sparray[lcol,i0+1])) then $
        newarray[rcol,i0]=newflux
      testlam=sparray[lcol,i1]
      if (abs(testlam-lambda) lt abs(testlam-sparray[lcol,i1-1])) then $
        newarray[rcol,i1]=newflux

    endelse

  endif

endfor

return,newarray
END
