FUNCTION spmax,inarray,l0,l1,lmax,RJ=rj,NUFNU=nufnu,FLAM=flam,ORDER=order

;  7 May 13 added an order keyword to limit the data to a single order
; 19 Jun 07 added lmax optional parameter as in spmin.pro
; 19 Oct 06 added /nufnu keyword for nu * F_nu (or lam * F_lam) units
; 22 Jun 06 added /rj keyword to compute max in RJ units
; 16 Mar 05 modified to set l0 and l1 if not in function call
; 31 Jan 05 created
;
; spmax returns the maximum flux between l0 and l1
;
; INPUT
;   inarray - input spectral array
;   l0,l1   - range to check for max, l1 must be greater than l0
;   lmax    - if given, used to return wavelength at maximum
; OUTPUT - returns max flux between l0 and l1

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
  maxflux=max(flux[index],maxspot)
  lmax=lambda[lcol,index[maxspot]]
  RETURN,maxflux
endif else begin
  print,'Error in spmax.  No data in range ',l0,' - ',l1
  RETURN,0
endelse

END
