PRO spfdivide,infile1,infile2,outfile

; 27 Jan 05 created
;
; INPUT
;   infile1 - spectral FITS file with numerator
;   infile2 - spectral FITS file with denominator
;   outfile - destination
; OUTPUT - writes ratioed spectrum in outfile
;
; no checking of array sizes

sp1=readfits(infile1,hdr)
sp2=readfits(infile2,hdr)

sp=spdivide(sp1,sp2)
writefits,outfile,sp,hdr

END
