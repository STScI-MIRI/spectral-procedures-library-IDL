FUNCTION spavg,sparray,l0,l1,LAM0=lam0,ERROR=error,ORDER=order,RJ=rj

; 19 Oct 05 - added lam0 keyword to return mean wavelength
;  9 Mar 05 - minor bug repair
; 29 Oct 04 - added /rj keyword
; 12 Oct 04 - troubleshooting; no changes
;  1 Apr 04 - added keyword error to return uncertainty in mean
;  4 Mar 04 - added keyword order to limit calculation to one order if needed
; 31 Dec 03 - fixed bug for overlapping segments
; 30 Nov 03 - created (based on the swsmake procedure segmentavg.pro)
;
; spavg finds the mean flux in the spectrum between l0 and l1
; the uncertainty in the mean is not propagated, but calculated from scratch
;
; INPUT
;   sparray - spectral data array
;   l0, l1  - starting and stopping wavelengths for averaging
;           - if l0 = l1 = 0.0, then average is over full data range
;   error   - optional parameter for return of error
;   order   - optional parameter to limit average to a specific order
;   rj      - optional keyword to force average in Rayleigh-Jeans flux units

lcol=0 & fcol=1 & ocol=3
if (keyword_set(order) eq 1) then begin
  idx=where(sparray[ocol,*] eq order)
  if (max(idx) gt -1) then begin
    l=reform(sparray[lcol,idx])
    if (keyword_set(rj) eq 0) then f=reform(sparray[fcol,idx]) $
      else f=reform(sparray[fcol,idx]*sparray[lcol,idx]^2)
  endif else begin
    print,'Error.  No data found for order ',order
    stop
  endelse
endif else begin
  l=reform(sparray[lcol,*])
  if (keyword_set(rj) eq 0) then f=reform(sparray[fcol,*]) $
    else f=reform(sparray[fcol,*]*sparray[lcol,*]^2)
endelse

len=long(n_elements(l))

if (l0 eq 0.0 and l1 eq 0.0) then begin ; if both=0, set to full range
  l0=min(l)
  l1=max(l)
endif ; shortened block 19 Oct 05
idx=where(l ge min([l0,l1]) and l le max([l0,l1]))

error=0.0
if (max(idx) gt -1) then begin
  num=float(n_elements(idx))
  avg=total(f[idx])/num
  lam0=mean(l[idx])
  if (num gt 1) then begin
    error = sqrt( total((f[idx]-avg)^2)/((num-1)*(num)) ) 
  endif
  return,avg
endif else return,0.0

END
