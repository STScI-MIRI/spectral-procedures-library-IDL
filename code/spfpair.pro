PRO spfpair,file_a,file_b,file_out,PLOT=plot,NOPAIR=nopair,_extra=e

; 19 Mar 09 now calling spplot with /error set
;  8 Apr 08 adding file name as title to plot
; 16 Jun 06 added NOPAIR keyword, if set, then finds mean instead
; 18 Dec 04 added PLOT keyword; default is now not to plot
;  9 Dec 04 created by GS
;
; wrapper for sppair.pro
; reads in two spectral FITS files, passes them to sppair.pro
; writes output as spectral FITS file using header from first input file

spa=readfits(file_a,hdr,/silent)
spb=readfits(file_b,/silent)

if (keyword_set(nopair) eq 0) then sp=sppair(spa,spb,_extra=e) $
  else sp=spmean(spa,spb)
if (keyword_set(plot) ne 0) then begin
  spplot,sp,title=file_out,/error
  zq=get_kbrd(1)
endif

writefits,file_out,sp,hdr

END
