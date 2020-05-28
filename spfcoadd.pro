PRO spfcoadd,infile,outfile,_extra=e

; 15 Dec 08 modified to pass input arguments to spcoadd
;  5 Mar 08 newly created
;           old spfcoadd is now spcoaddbat
;
; spfcoadd loads a file with FITS file names, calls spcoadd
;   if /plot set, then iteratively moves through the file list
;     calling spcoadd with one more file each time
;   if outfile specified, then writes coadded spectrum to a spectral FITS file
;
; INPUT
;   infile  - file with one spectral FITS file name per line
;   outfile - if specified, then coaded spectrum will be written to outfile

openr,fi,infile,/get_lun
line=infile

; open infile and load files within it into a single string array

count=0
while (not eof(fi)) do begin
  readf,fi,line
  if (count eq 0) then list=line else list=[list,line]
  count=count+1
endwhile

; if plotting, then plot first spectrum and make repeated calls to spcoadd, 
; adding one spectrum at a time and replotting

if (keyword_set(plot) eq 1) then $
  spplot,readfits(list[0],/silent),tit=list[0],_extra=e

sp=spcoadd(list,hdr,_extra=e)

; write output FITS file if requested

if (n_elements(outfile) ne 0) then writefits,outfile,sp,hdr
    
close,/all
END
