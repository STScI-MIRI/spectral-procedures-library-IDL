FUNCTION spse,sp,csp,LAM=lam,DLAM=dlam,PLOT=plot,QUIET=quiet,_extra=e

;  3 Jul 08, now scaling csp so that it can be returned to calling procedure
; 13 Sep 07, added DLAM keyword (for measuring DEC), changed order of output
;            switched from use of spsum.pro to spint.pro
; 23 Dec 06, added /quiet keyword, also passing _extra keywords to spplot
;  3 Apr 05, added LAM as optional keyword (to input continuum-fitting range)
; 18 Nov 04, modified lam range from 7.6-8.8 to 7.6-8.6
; 23 Aug 04, as spfse, but receives spectral arrays, not files
; 23 Aug 04, renamed from se_class.pro now requires continuum file as input
; 12 Apr 04, created by M.J. Russell from C code by G.C. Sloan
;
;  takes an LRS spfits image, fits and subtracts aa continuum,
;  computes averages, at 10,11,12 microns
;  finds flux ratios 10/11, 10/12, corrected ratio 11/12 and contrast
;
;  INPUT
;    sp - a spectral data array for science spectrum
;    csp - ditto for the continuum spectrum, to be scaled to sp
;    lam - optional 2-element vector, gives range for continuum fitting
;          default defined below
;    dlam - optional 2-element vector, gives range to calculate DEC
;           default defined below
;    plot - optional keyword to plot sp, fitted csp, and resulting dust
;    quiet - turns off prompts for user to continue when plotting
;  OUTPUT
;    returns a vector with following elements:
;      DEC & error F10/F11 & error F10/F12 & error, 'F11/F12' & error
;    also returns scaled version of continuum in csp variable

if (not keyword_set(plot)) then plot=0

; check lam and dlam keywords, set to default if not set

if (n_elements(lam)  ne 2) then  lam=[7.60,8.60] 
if (n_elements(dlam) ne 2) then dlam=[7.50,14.50]

; normalize continuum spectrum to science spectrum
;   over the range lam[0:1] or 7.60-8.80 um (7.67-8.74 um - first 7 pixels)
;   and subtract it

factor = spavg(sp,lam[0],lam[1])/spavg(csp,lam[0],lam[1])
csp = sptimes(csp,factor)
dust=spadd(sp,csp,/minus)

if (plot eq 1) then begin
   spplot,sp,th=1,_extra=e
   spplot,csp,th=2,li=1,/over,_extra=e
   spplot,dust,th=2,/over,_extra=e
   if (not keyword_set(quiet)) then print,'Hit a key to continue.'
   zchar=get_kbrd(1)
endif

; calculate fluxes at 10,11,12 um and their errors
; find ratios and propagate errors

f10 = spavg(dust,9.756,10.185,error=err10)
f11 = spavg(dust,10.859,11.242,error=err11)
f12 = spavg(dust,11.850,12.199,error=err12)

f10_11 = f10/f11
f10_12 = f10/f12
err10_11 = sqrt( (f10_11)^2*( (err10/f10)^2 + (err11/f11)^2 ) )
err10_12 = sqrt( (f10_12)^2*( (err10/f10)^2 + (err12/f12)^2 ) )

; move ratio back to power law, calculate 11/12 ratio and error

corratio = mv2curve(f10_11,f10_12)
f11_12 = corratio[1]/corratio[0]
err11_12 = sqrt( f11_12^2*( (err10_12/f10_12)^2 + (err10_11/f10_11)^2 ) )

;calculate dust contrast

d_err=1
s_err=1
;d_sum = spsum(dust, dlam[0],dlam[1],error=d_err)
;s_sum = spsum(csp,dlam[0],dlam[1],error=s_err)
;dcon = d_sum/s_sum
;decerr = sqrt( dcon^2*( (d_err/d_sum)^2 + (s_err/s_sum)^2 ) )

d_sum = spint(dust, dlam[0],dlam[1],/error)
s_sum = spint(csp,dlam[0],dlam[1],/error)
dcon = d_sum[0]/s_sum[0]
decerr = sqrt( dcon^2*( (d_sum[1]/d_sum[0])^2 + (s_sum[1]/s_sum[0])^2 ) )

; return dust em constrast, error, f10/f11, err., f10/f12, err., f11/f12, err.

return,[dcon,decerr,f10_11,err10_11,f10_12,err10_12,f11_12,err11_12]

end
