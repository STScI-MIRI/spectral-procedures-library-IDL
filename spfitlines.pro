FUNCTION spfitlines,sparray,lambda,degree,SUBTRACT=subtract

; 30 Mar 05 created
;
; spfitlines fits a series of line segments to selected regions of sparray

; INPUT
;   sparray - input spectrum, must have wavelength and flux (col 0 and 1)
;   lambda  - an array of wavelengths, must be divisible by four
;   subtract - optional keyword to subtract the fitted polynomial 
;             from the spectrum
; OUTPUT - returns a spectral data array with the line segments
;          unless the subtract keyword is set, in which case the
;          difference between the input spectrum and the polynomial is returned

lcol=0 & fcol=1

; check input lambda array, make sure it is divisible by four
; then load l0 and l1

nlam=n_elements(lambda)
nset=nlam/4
if (nset*4 ne nlam) then begin
  print,'Warning in spfitlines.  Input wavelengths must be in sets of four.'
endif

; load l0 and l1, starting and stopping wavelength ranges for polynomial fitting

l0=fltarr(nset) & l1=l0 & l2=l0 & l3=l0
for i=0,nset-1 do begin
  l0[i]=lambda[4*i]
  l1[i]=lambda[4*i+1]
  l2[i]=lambda[4*i+2]
  l3[i]=lambda[4*i+3]
endfor

; find indices within ranges l0[0]:l1[0], l2[0]:l3[0], l0[1]:l1[1], etc.

deg=2
continuum=sparray

for i=0,nset-1 do begin

; fit a line to the data l0:l1, l2:l3
; use it to replace data 

  fit_idx=where( (sparray[lcol,*] ge l0[i] and sparray[lcol,*] le l1[i]) $
             or  (sparray[lcol,*] ge l2[i] and sparray[lcol,*] le l3[i]) )
  mod_idx=where(  sparray[lcol,*] ge l0[i] and sparray[lcol,*] le l3[i]  )
  if (max(fit_idx) gt -1 and max(mod_idx) gt -1) then begin
    scc=svdfit(reform(sparray[lcol,fit_idx]),reform(sparray[fcol,fit_idx]),deg)
    continuum[fcol,mod_idx]=poly(reform(continuum[lcol,mod_idx]),scc)
  endif else begin
    print,'Warning in spfitlines.  Check wavelengths.  No data in range.'
  endelse
endfor

if (n_elements(subtract) eq 0) then begin
  RETURN,continuum
endif else begin
  RETURN,spadd(sparray,continuum,/minus)
endelse

END
