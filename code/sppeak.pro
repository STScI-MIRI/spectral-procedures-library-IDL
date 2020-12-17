FUNCTION sppeak,sp,l0,l1

; 25 Jan 05 created
; 
; given an input spectrum, report wavelength of maximum flux
; if l0 and l1 specified, limit search for max to that wavelength range

lcol=0 & fcol=1 ; ocol not used (yet)
l=sp[lcol,*] & f=sp[fcol,*]

if (n_elements(l0) gt 0 and n_elements(l1) gt 0) then begin
  idx=where(l ge l0 and l le l1)
  if (max(idx) lt 0) then begin
    print,'Warning in sppeak.pro - no data in range ',l0,l1,'.  Returning 0.'
  endif
  ll=l[idx] & ff=f[idx]
endif else begin
  ll=l      & ff=f
endelse

fmax=max(ff,idx_max)

RETURN,ll[idx_max]
END
