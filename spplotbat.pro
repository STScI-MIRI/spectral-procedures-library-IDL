PRO spplotbat,inlist,xrange=xrange,yrange=yrange,lambda=lambda,spline=spline,same=same,contfit=contfit,skip=skip,vline=vline,irs=irs,_extra=e

; 20 Jul 14 modifying to call rd_spfits or rd_sptbl for .fits and .tbl suffixes
; 27 Dec 13 skipping lines with zero length
; 24 Jun 13 added irs keyword to force readfits call (for most IRS data)
; 11 Nov 12 updating contfit and contflag (differs from conflag)
; 22 Oct 12 replaced readfits call with rd_spfits call (for old spec FITS files)
; 26 Aug 12 added vline keyword to overplot vertical lines at listed wavelengths
;  2 Aug 12 added the skip keyword to skip header lines in the list file
; 19 Jul 12 added spline and contfit keywords
;  2 Mar 12 added same keyword to turn off autoscaling of yrange after 1st plot
; 26 Oct 11 added the lambda keyword, which invokes a continuum subtraction
; 22 Jan 10 created from scratch
;
; accepts either
; (1) an array of strings containing a list of spectral FITS files
; (2) a file containing the list of spectral FITS files
; and displays them, using keyboard input and dobat.pro to navigate the list
;
; lam keyword - sets wavelengths for continuum subtraction with spcont.pro
; 
; note that yrange is implemented poorly at this point - the user needs to
; toggle up and down in the input list to reset the range correctly
;
; INPUT
;   xrange  - passed to spplot
;   yrange  - passed to spplot
;   lambda  - only four values allowed - will call spcont and display diff
;   spline  - if set, the displays spline-subtracted spectra (overrides lambda)
;   contfit - only used if lambda or spline used - plots continuum fit, not diff
;   same    - prevents autoscaling of y range from spectrum to spectrum
;   skip    - tells spplotbat to skip header lines in the list file

istop=n_elements(inlist)-1

; check keywords
; set conflag (whether removing a continuum or not)
; and check lambda, spline, vline keywords

conflag=0
if (keyword_set(lambda) ne 0) then begin 
  if (n_elements(lambda) eq 4) then conflag=1 else $
    print,'Warning.  Must give FOUR elements for lambda.'
endif
if (keyword_set(spline) ne 0) then begin
  if (n_elements(spline) gt 2) then conflag=2 else $
    print,'Warning.  Must give at least THREE elements for spline.'
endif
if (keyword_set(contfit) eq 0) then contflag=0 else begin
  contflag=1
  if (conflag eq 0 and contfit eq 1) then contflag=0
endelse

nlines=n_elements(vline)

; load the input list

if (istop eq 0) then begin ; read in the input file list, load to string array
  openr,fin,inlist,/get_lun
  count=0 & instr=' '
  if (keyword_set(skip) ne 0) then for i=0,skip-1 do readf,fin,instr
  while (not eof(fin)) do begin
    readf,fin,instr
    if (strlen(instr) gt 0) then begin
      split=strsplit(instr,' ',/extract)
      if (count eq 0) then file=split[0] else file=[file,split[0]]
      count=count+1
    endif
  endwhile
endif else file=inlist ; inlist is already a string array
count=count-1 ; for future use as a stop

; handle xrange and yrange input specially
; everything else just passed to spplot calls (or not)

if (n_elements(xrange) ne 2) then xrange=0

case n_elements(yrange) of
  2 : range=yrange
  else : range=0
endcase

goflag=1
counter=0
print,"Type a command after the spectrum is plotted - 'h' for help."

while (goflag eq 1) do begin

; display image (if it can be read)

;  if (keyword_set(irs) eq 0) then 
   split=strsplit(file[counter],'.',/extract)
   case split[n_elements(split)-1] of
     'tbl' : data=rd_sptbl(file[counter])
     else  : data=rd_spfits(file[counter],/silent) ; else $
;    data=readfits(file[counter],/silent)
   endcase

  if (n_elements(data) eq 1) then begin
    if (data eq -1) then print,'Error.  Cannot read file ',file[counter],' .'
    if (data ne -1) then print,'Error.  Illegal file.'
    print,'No display'
  endif else begin

;   if conflag set, subtract continuum from data

    if (conflag ne 0) then begin
      case conflag of
        1 : cont=spcont(data,lambda[0],lambda[1],lambda[2],lambda[3],_extra=e)
        2 : cont=spspline(data,spline,_extra=e)
      endcase
      if (contflag eq 0) then data=spadd(data,cont,/minus)
    endif

    spplot,data,xrange=xrange,yrange=range,title=file[counter],_extra=e
    if (contflag ne 0) then spplot,cont,/over,li=1
    for i=0,nlines-1 do oplot,[vline[i],vline[i]],range,li=1
  endelse

; pass input to dobat, update counter, goflag, range, zoom, or xrange, continue

  oldcounter=counter
  dobat,get_kbrd(1),counter,count,$
    goflag=goflag,range=range,zoom=zoom,xrange=xrange

; if counter has changed and same keyword not set, then
;   reset yrange so that next plot is autoscaled

  if (counter ne oldcounter and keyword_set(same) eq 0) then range=0 

endwhile

free_lun,fin
END
