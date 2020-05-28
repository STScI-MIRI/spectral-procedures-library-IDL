FUNCTION spchisq,sparray,l0,l1,ORDER=order

; 26 Apr 06 if l0, l1 undefined, then set to min and max of wavelength range
; 29 Apr 04 created by modifying spavg.pro
;
; spchisq calculates the square of the sum of sparray between l0 and l1
; in general, one should pass spchisq a difference 
; spavg finds the mean flux in the spectrum between l0 and l1
; the uncertainty in the mean is not propagated, but calculated from scratch
; if l0 = l1 = 0.0, then average is over full data range

lcol=0 & fcol=1

if (keyword_set(order) eq 1) then begin
  index=where(sparray[3,*] eq order)
  if (max(index) gt -1) then begin
    l=reform(sparray[lcol,index])
    f=reform(sparray[fcol,index])
  endif else begin
    print,'Error.  No data found for order ',order
    stop
  endelse
endif else begin
  l=reform(sparray[lcol,*])
  f=reform(sparray[fcol,*])
endelse

len=long(n_elements(l))

; if both l0 and l1 not defined or zero, then set l0-l1 to full wavelength range

if (keyword_set(l0) eq 0 and keyword_set(l1) eq 0) then begin 
  l0=min(l)
  l1=max(l)
endif 

; if l0>l1, flip l0 and l1

if (l0 gt l1) then begin     
  ldum=l0
  l0=l1
  l1=ldum
endif

idx=where(l ge l0 and l le l1)

sumsq=0.0
if (max(idx) gt -1) then begin
  for i=0,n_elements(idx)-1 do sumsq = sumsq +f[idx[i]]^2.0
endif

return,sumsq

END
