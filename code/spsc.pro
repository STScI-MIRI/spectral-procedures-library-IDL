FUNCTION spsc,sp,csp,LAM=lam,PLOT=plot

;  3 Apr 05 added LAM keyword to change default wavelengths for continuum fit
; 21 Jan 05 modified to assume csp=0 if not passed
; 18 Nov 04 modified lam range from 7.6-8.8 to 7.6-8.6
;  5 Aug 04 printing out dust em. contrast now
; 24 Aug 04 renamed from spse, modified lambda ranges, no mv2curve
; 23 Aug 04 as spfse, but receives spectral arrays, not files
; 23 Aug 04 renamed from se_class.pro now requires continuum file as input
; 12 Apr 04 created by M.J. Russell from C code by G.C. Sloan
;
;  takes an LRS spfits image, fits and subtracts aa continuum,
;  computes average flux at 9.7, 10.3, 10.7 (F_a, F_b, F_c)
;  finds flux ratios F_a/F_b, F_a/F_c, dust emission contrast (DEC)
;  DEC = ratio of dust/continuum from 7.5 to 14.5 um
;  does not find corrected ratio F_a/F_c
;
; INPUT
;   sp - spectral data array
;   csp - matching spectral data array with continuum spectrum
;       - if csp not provided, continuum assumed to be zero
;   lam - optional keyword to change default wavelengths to fit continuum
; OUTPUT - returns data and errors for F_a/F_b, F_a/F_c, DEC
;          F_a, F_b, F_c = equivalent flux in three windows
;          DEC = dust emission contrast

if (not keyword_set(plot)) then plot=0

; normalize continuum spectrum to science spectrum
;   over the range 7.60-8.80 um (7.67-8.74 um - first 7 pixels)
;   and subtract it

if (n_elements(lam) eq 0) then begin
  l0=7.60 & l1=8.60
endif else begin
  l0=lam[0] & l1=lam[1]
endelse

if (n_elements(csp) gt 0) then begin
  factor = spavg(sp,l0,l1)/spavg(csp,l0,l1)
  fitsp = sptimes(csp,factor)
  dust=spadd(sp,fitsp,/minus)
endif else begin
  fitsp=sp
  dust=sp
endelse

if (plot eq 1) then begin
   !x.range=0 & !x.style=0 & !y.range=0 & !y.style=0
   spplot,sp,th=1
   spplot,fitsp,th=2,li=1,/over
   spplot,dust,th=2,/over
   print,'Hit a key to continue.'
   zchar=get_kbrd(1)
endif

; calculate fluxes at 10,11,12 um and their errors
; find ratios and propagate errors

f_a = spavg(dust,9.50,9.90,error=err_a)
f_b = spavg(dust,10.05,10.40,error=err_b)
f_c = spavg(dust,10.45,10.85,error=err_c)

f_a_b = f_a/f_b
f_a_c = f_a/f_c
err_a_b = sqrt( (f_a_b)^2*( (err_a/f_a)^2 + (err_b/f_b)^2 ) )
err_a_c = sqrt( (f_a_c)^2*( (err_a/f_a)^2 + (err_c/f_c)^2 ) )

;calculate dust contrast
                                                                                
d_err=1
s_err=1
d_sum = spsum(dust, 7.5,14.5,error=d_err)
s_sum = spsum(fitsp,7.5,14.5,error=s_err)
dcon = d_sum/s_sum
decerr = sqrt( dcon^2*( (d_err/d_sum)^2 + (s_err/s_sum)^2 ) )

; return f_a/f_b, error, f_a/f_c, error 

return,[f_a_b,err_a_b,f_a_c,err_a_c,dcon,decerr]

end
