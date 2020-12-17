FUNCTION spsum,sp,l0,l1,error=ERROR

; 13 Sep 07 tried to debug after finding discrepancies with similar spint.pro
;           repaired bug with indexing in for loop
; 18 Oct 05 sort spectral array, general delta lambda before integrating
; 12 Apr 04 modifications to comments by GS
; 11 Apr 04 created by MJ Russell
;
; numerically integrates a spectrum from l0 to l1
; also returns rms error if passed

sparray=spsort(sp,/nodupe)
l=reform(sparray[0,*])
f=reform(sparray[1,*])
e=reform(sparray[2,*])

; convert to F-lambda

f = 3.0e-12*f/(l^2.0)
e = 3.0e-12*e/(l^2.0)

; generate delta Lam

dellam=l-shift(l,1)
dellam[0]=dellam[1]

; set index of wavelength elements to integrate

idx=where(l ge min([l0,l1]) and l le max([l0,l1]))

; integrate and return flux

ftotal=0.0
etotal=0.0
error=0.0
if (max(idx) gt -1) then begin
  for i=0,n_elements(idx)-2 do begin
    ftotal = ftotal + 0.5*(f[idx[i+1]] + f[idx[i]])*dellam[idx[i]]
    etotal = etotal + (e[idx[i]])^2 
  endfor
  error = sqrt(etotal/(n_elements(idx)-1))
endif

return,ftotal

END
