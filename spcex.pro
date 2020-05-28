FUNCTION spcex,insp,l0,l1,l2,l3,lamc,$
  DEGREE=degree,MODE=mode,FLUX=flux,RATIO=ratio,ORDER=order,PLOT=plot,_extra=e

;  7 May 13 added order keyword
; 11 Feb 11 now passing _extra to spcont call in case the /RJ flag is set
; 27 Feb 07 moved all error propagation into spex.pro
; 12 Sep 05 now calling spex using mode keyword
; 31 Aug 05 renamed as spcex, corrected errors for /flux and /ratio calls
; 29 Aug 05 added lamc argument, made degree an optional parameter
; 26 Mar 05 returned uncertainty now includes uncertainty in continuum
; 16 Mar 05 improved plotting with _extra=e
;  5 Mar 05 created as spcexw
;
; calls spcont to fit a continuum between l0 and l3
; then extracts equivalent width between l1 and l2

; INPUT
;   l0,l1   wavelength range on blue side of feature
;   l2,l3   wavelength range on red side of feature
;   lamc    wavelength 
;   degree  keyword, degree of polynomial; defaults to 2 (line)
;   mode    keyword, if set to 1, calls spex with /flux
;                    if set to 2, calls spex with /ratio
;   flux    keyword to explicitly set /flux in spex call
;   ratio   keyword to explicitly set /ratio in spex call
;   order   keyword to limit analysis to a single order of the input spectrum

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check keywords

if (keyword_set(degree) eq 0) then degree=2
if (keyword_set(mode) eq 0) then begin
  mode=0
  if (keyword_set(flux) eq 1) then mode=1
  if (keyword_set(ratio) eq 1) then mode=2
endif

if (keyword_set(order) ne 0) then begin
  idx=where(insp[ocol,*] eq order)
  if (idx[0] gt -1) then sp=insp[*,idx] else sp=insp
endif else sp=insp

; determine continuum spectrum, plot if /plot set

spc=spcont(sp,l0,l1,l2,l3,degree,sigratio,_extra=e)

if (keyword_set(plot) eq 1) then begin
  spplot,sp,_extra=e
  spplot,spc,_extra=e,/over
  ymin=spmin(sp) & ymax=spmax(sp)
  if (ymin gt 0) then ymin=0.5*ymin else ymin=2.0*ymin
  if (ymax gt 0) then ymax=2.0*ymax else ymax=0.5*ymax
  oplot,[l0,l0],[ymin,ymax],_extra=e
  oplot,[l1,l1],[ymin,ymax],_extra=e
  oplot,[l2,l2],[ymin,ymax],_extra=e
  oplot,[l3,l3],[ymin,ymax],_extra=e
  zq=get_kbrd(1)
endif

; call spex to find equivalent width, flux, or flux ratio

eqw=spex(sp,spc,l1,l2,lamc,sigratio,mode=mode,_extra=e)

; all error propagation is now done in spex

return,eqw
END
