PRO speplot,in_data,ORDER=order,_extra=e,rj=rj,flam=flam,lamflam=lamflam,nufnu=nufnu

; 12 Jul 10 added keyword flam for f_lam units (F_nu*c*lam^-2)
;           and keywords lam, lamflam, nufnu for lam f_lam units
;  7 Mar 08 now can read a FITS file
; 18 Jan 08 updated for loop counter to long for BIG spectra
; 28 Sep 07 added order keyword (!)
; 18 May 04 added keyword rj for Rayleigh-Jeans plot (F_nu*lam^2)
; 28 Nov 03 created
;
; plots error bars on a spectral plot
; expects col 0 = wavelength, 1 = flux, 2 = uncertainty, 3 = order

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check in_data to see if it's an array or a file name, if the latter, read it

if (n_elements(in_data) gt 1) then in_array=in_data else $
  in_array=readfits(in_data,/silent)

; check data and keywords

if (n_elements(in_array[*,0]) lt 3) then begin
  print,"Error.  Array must have 3 or more columns"
  stop
endif

if (keyword_set(order) eq 0) then data=in_array else begin
  idx=where(in_array[ocol,*] eq order)
  if (max(idx) gt -1) then data=in_array[*,idx] else begin
    print,'Error in speplot.  No data of order ',order
    stop
  endelse
endelse

if (keyword_set(rj)      eq 1) then   rjflag=1 else   rjflag=0
if (keyword_set(flam)    eq 1) then flamflag=1 else flamflag=0
if (keyword_set(lamflam) eq 1 or keyword_set(nufnu)   eq 1) $
                               then lamflag=1 else  lamflag=0

; plot

len=n_elements(data[lcol,*])
for i=long(0),len-1 do begin
  xx=[data[lcol,i],data[lcol,i]]
  yy=[data[fcol,i]-data[ecol,i],data[fcol,i]+data[ecol,i]]
  if (rjflag eq 1) then yy=yy*xx^2 $
  else if (flamflag eq 1) then yy=yy*3e-12/(xx^2) $
  else if (lamflag eq 1) then yy=yy*3e-12/xx
  oplot,xx,yy,_extra=e
endfor

END
