FUNCTION spmin,inarray,l0,l1,lmin,RJ=rj,NUFNU=nufnu,FLAM=flam,ORDER=order

;  7 May 13 updating to include changes to smpax.pro
; 19 Jun 07 added /nufnu and /rj keywords as in spmax.pro
;  3 Jan 06 added lmin optional parameter
; 16 Mar 05 created by modifying spmax.pro
;
; spmin returns the minimum flux between l0 and l1
;
; INPUT
;   inarray - input spectral array
;   l0,l1   - range to check for max, l1 must be greater than l0
;   lmin    - if given, used to return wavelength at minimum
; OUTPUT - returns min flux between l0 and l1

lcol=0 & fcol=1 & ocol=2

if (keyword_set(order) ne 0) then begin
  idx=where(inarray[ocol,*] eq order)
  if (idx[0] gt -1) then sparray=inarray[*,idx] else sparray=inarray
endif else sparray=inarray

lambda=sparray[lcol,*]
flux=sparray[fcol,*]

if (keyword_set(rj) eq 1) then flux=sparray[fcol,*]*lambda^2
if (keyword_set(flam) eq 1) then flux=sparray[fcol,*]*3e-12/lambda^2
if (keyword_set(nufnu) eq 1) then flux=sparray[fcol,*]*3e-12/lambda

if (n_elements(l0) eq 0) then l0=min(lambda)
if (n_elements(l1) eq 0) then l1=max(lambda)
index=where(lambda ge l0 and lambda le l1)

if (max(index) gt -1) then begin
  minflux=min(flux[index],minspot) 
  lmin=lambda[lcol,index[minspot]]
  RETURN,minflux
endif else begin
  print,'Error in spmin.  No data in range ',l0,' - ',l1
  RETURN,0
endelse

END
