FUNCTION spcont,sp,l0,l1,l2,l3,degree,sigratio,RJ=rj,PLOT=plot,_extra=e

; 23 Jul 12 addec plot keyword
; 11 Feb 11 added RJ mode to perform calculations in RJ units
; 26 Mar 05 added simple error calculation, returned as sigratio
;  5 Mar 05 created
;
; spcont fits a continuum to a spectrum
; outside of l0 and l3, the continuum matches the spectrum
; a polynomial is fit to the data in the ranges l0-l1 and l2-l3
; and applied from l0 to l3
;
; INPUT
;  l0,l1     wavelength range on blue side of feature
;  l2,l3     wavelength range on red side of feature
;  degree    degree of polynomial; defaults to 2 (line)
;  sigratio  optional parameter to return uncertainty in continuum
;            (as a fraction of continuum flux) between l1 and l2

; prepatory definitions and checks

lcol=0 & fcol=1 & ocol=3
if (keyword_set(degree) lt 1) then degree=2
if (keyword_set(rj) eq 0) then spreturn=sp else spreturn=sprj(sp)

; load wavelength array and set stops

l=reform(sp[lcol,*])
idx=where( (l ge l0 and l le l1) or (l ge l2 and l le l3) )
idxfull=where( l ge l0 and l le l3)

if (max(idx) gt -1) then begin
  lam=reform(spreturn[lcol,idx])
  flux=reform(spreturn[fcol,idx])
  cc=svdfit(lam,flux,degree)
  spreturn[fcol,idxfull] = poly(reform(spreturn[lcol,idxfull]),cc)

; calculate error in continuum flux (uncertainty in mean)

   fsig = stddev(spreturn[fcol,idx]-sp[fcol,idx]) / sqrt(n_elements(idx))

; scale error to average continuum flux

  sigratio = fsig / ( total(sp[fcol,idx]) / n_elements(idx) )

endif else begin
  print,'Warning.  No data in fitting range, returning input spectrum'
  sigratio=0.0
endelse

if (keyword_set(plot) ne 0) then begin
  spplot,sp,_extra=e
  spplot,spreturn,/over,li=1,_extra=e
  zq=get_kbrd(1)
endif

if (keyword_set(rj) ne 0) then spreturn=sprjout(spreturn)

return,spreturn
END
