PRO spfnorm,infile,outfile,norms

; 29 Dec 04 created
;
; wrapper for spnorm, see it for comments
;
; INPUT
;   infile  - input spectral FITS file
;   outfile - output spectral FITS file
;   norms   - array of scalar multiplicative corrections
;             indices of corrections must match order numbers
;             this usually requires leading zeroes (up to 11!)
; OUTPUT - writes results to outfile

spin=readfits(infile,hdr,/silent)
spout=spnorm(spin,norms)
writefits,outfile,spout,hdr

END
