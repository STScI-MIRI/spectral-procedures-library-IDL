PRO spfadd,infile1,infile2,outfile,_extra=e

; 16 Jun 06 created
;
; INPUT
;   infile1 - spectral FITS file 
;   infile2 - spectral FITS file 
;   outfile - destination file for sum or difference
; OUTPUT - writes summed or differneced spectrum in outfile
;
; no checking of array sizes

sp1=readfits(infile1,hdr)
sp2=readfits(infile2,hdr)

sp=spadd(sp1,sp2,_extra=e)
writefits,outfile,sp,hdr

END
