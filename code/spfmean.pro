PRO spfmean,file1,file2,outfile,PLOT=plot,_extra=e

; 17 Sep 05 created by modifying spfjoin
;
; spfmean averages two spectral FITS files
; they must have the same number of columns
; the header from the first file will be used as the output FITS header
; assumed that col 3 = order number

sp1=readfits(file1,hdr,/silent)
sp2=readfits(file2,/silent)

sp=spmean(sp1,sp2,_extra=e)

wr_spfits,outfile,sp,-1,fits_hdr=hdr

if (keyword_set(plot) eq 1) then begin
  spplot,sp,_extra=e
  speplot,sp,_extra=e
  zq=get_kbrd(1)
endif

END
