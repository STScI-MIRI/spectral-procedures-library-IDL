FUNCTION spcoadd,filelist,hdr,spike=SPIKE,normal=NORMAL,plot=PLOT,_extra=e

; 30 Mar 20 header now correctly describes procedure (doesn't write to disk)
; 15 Dec 08 added normal keyword to normalize all data before coadding
;  5 Mar 08 modified to a function, returns data, writes nothing to disk now
; 26 Dec 03 add spike keyword
; 22 Dec 03 renamed to spcoadd.pro and moved to my procs directory
; 17 Dec 03 created as coadd.pro

; reads a list of spectral files, coadds them, and returns coadded result
;  as a spectral FITS file (but with no NSEG header keywords)
; assumes all are on the same grid (for now)
;
; INPUT
;   filelist - string array with list of spectral FITS file to read
; OUTPUT
;   hdr      - if passed, used to return FITS header of first file read
;   returns coadded spectral data array

lcol=0 & fcol=1 & ecol=2 & ocol=3

if (keyword_set(spike) eq 0) then spikeflag=0 else spikeflag=1
if (keyword_set(plot) eq 0) then plotflag=0 else plotflag=1

print,keyword_set(normal)

nfiles=n_elements(filelist)

for i=0,nfiles-1 do begin

; read in file names and load these one by one with readfits
                                                                                
  sp=readfits(filelist[i],tmphdr,/silent)
  if (keyword_set(normal) ne 0) then sp=spdivide(sp,spmax(sp)) ; normalize 
  if (i eq 0) then hdr=tmphdr

; plot

  if (plotflag eq 1) then spplot,sp,_extra=e

; load data arrays

  lam=reform(sp[lcol,*])
  flux=reform(sp[fcol,*])
  order=reform(sp[ocol,*])
  nlen=n_elements(lam)

  if (spikeflag eq 1) then newflux=fix_spikes(flux,range=4,threshold=2.0) $
    else newflux=flux
  if (i eq 0) then fluxarray=newflux else fluxarray=[[fluxarray],[newflux]]

endfor

; determine median flux and uncertainty in mean

flux=fltarr(nlen)
unc=fltarr(nlen)

for j=0,nlen-1 do begin
  flux[j]=median(fluxarray[j,*])
  unc[j]=stddev(fluxarray[j,*])/sqrt(nfiles)
endfor

; load output array

outarray=fltarr(4,nlen)
outarray[lcol,*]=lam
outarray[fcol,*]=flux
outarray[ecol,*]=unc
outarray[ocol,*]=order

return,outarray
END
