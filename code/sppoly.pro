FUNCTION sppoly,sparray,lambda,degree,SUBTRACT=subtract

; 31 Jan 05 created
;
; sppoly fits a polynomial to the selected regions of sparray

; INPUT
;   sparray - input spectrum, must have wavelength and flux (col 0 and 1)
;   lambda  - an array of wavelengths, must have an even number
;   degree  - the degree of the polynomial to fit to sparray
;   subtract - optional keyword to subtract the fitted polynomial 
;             from the spectrum
; OUTPUT - returns a spectral data array with the fitted polynomial
;          unless the subtract keyword is set, in which case the
;          difference between the input spectrum and the polynomial is returned

lcol=0 & fcol=1

; check input lambda array, make sure it is even, then load l0 and l1

nlam=n_elements(lambda)
npair=nlam/2
if (npair*2 ne nlam) then begin
  print,'Warning in sppoly.  Must give an even number of wavelengths.'
endif

; load l0 and l1, starting and stopping wavelength ranges for polynomial fitting

l0=fltarr(npair) & l1=l0
for i=0,npair-1 do begin
  l0[i]=min([lambda[2*i],lambda[2*i+1]])
  l1[i]=max([lambda[2*i],lambda[2*i+1]])
endfor

; find indices within ranges l0[0]:l1[0], l0[1]:l1[1], etc.

paircount=0
for i=0,npair-1 do begin
  pair_idx=where(sparray[lcol,*] ge l0[i] and sparray[lcol,*] le l1[i])
  if (max(pair_idx) gt -1) then begin
    if (paircount eq 0) then idx=pair_idx else idx=[idx,pair_idx]
    paircount=paircount+1
  endif
endfor

; find polynomial and calculate continuum array

continuum=sparray
if (max(idx) gt -1) then begin
  scc=svdfit(reform(sparray[lcol,idx]),reform(sparray[fcol,idx]),degree)
  continuum[fcol,*]=poly(reform(continuum[lcol,*]),scc)
endif else print,'Error in sppoly.  No data in specified ranges.'

if (n_elements(subtract) eq 0) then begin
  RETURN,continuum
endif else begin
  RETURN,spadd(sparray,continuum,/minus)
endelse

END
