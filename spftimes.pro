PRO spftimes,infile,scale,outfile,_extra=e

;  3 Jan 06 added _extra=e
; 30 Oct 04 created
;
; INPUT
;   infile - input spectral FITS file
;   scale  - scalar quantity to multiply the spectrum by
;   outfile - destination
; OUTPUT - writes scaled spectrum in outfile

sp=readfits(infile,hdr)
sp=sptimes(sp,scale,_extra=e)
writefits,outfile,sp,hdr

END
