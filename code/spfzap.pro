PRO spfzap,infile,outfile,lam1,lam2,_extra=e

; 22 Apr 10 created
;
; INPUT
;   infile  - spectral FITS file to be trimmed (zapped)
;   outfile - destination for trimmed file
;   lam1    - starting wavelength of interval to be removed
;   lam2    - ending wavelength of interval to be removed
; OUTPUT - zaps requested wavelength and writes result to outfile

spin=readfits(infile,hdr,/silent)
spdata=spzap(spin,lam1,lam2,_extra=e)
wr_spfits,outfile,spdata,-1,FITS_HDR=hdr

END
