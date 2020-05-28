PRO spfjoin,file1,file2,outfile,BUMP=bump,BUMP2=bump2,NOBUMP=nobump,CROSS=cross,PLOT=plot,_extra=e

; 29 Jul 14 adding columns as needed if files differ
;  3 Mar 14 added bump2 keyword to fix the bump in order numbers for 2d spectrum
;  9 Dec 11 again modifying how orders numbers are bumped
;           del is now max of nseg1, nseg2, max1-min2+1
;  6 Dec 11 added cross keyword to call crossradec to improve coordinates
; 30 Nov 11 added bump and nobump keywords, now accounting for missing files
;           returning to older determination of order increment
;           was max1-min2+1 - not sure why I switched
; 17 Sep 05 order increment = max number of segments in input spectra
; 14 Mar 05 added _extra keyword to pass to spplot call
; 29 Dec 04 added a get_kbrd request after the plot, include file name in plot
; 18 Dec 04 added keyword PLOT to turn off plotting; default is now not to plot
;  8 Mar 04 renamed from spjoin to spfjoin
;  2 Feb 04 pass header to wr_spfits, overlooked this before
; 27 Dec 03 created as spjoin
;
; spfjoin combines as sequential segments two spectral FITS files
; they must have the same number of columns
; the header from the first file will be used as the output FITS header
; the order numbers in file 2 are incremented unless nobump is explicitly set
;   this is true whether or not file 1 is valid
; the order numbers in file 1 are incremented only:
;   if file 2 is not specified and bump keyword set
;   note that file 1 can never be incremented if file 2 is specified
; assumed that col 3 = order number - should be read from FITS header
;
; INPUT
;   file1   - first FITS file to be joined
;   file2   - second FITS file to be joined
;   bump    - optional keyword to bump order numbers if only one file read 
;   bump2   - optional keyword to fix bump in order number for file2
;   nobump  - optional keyword to NOT bump orders on file2 - overrides bump
;   plot    - optional keyword to plot joined spectrum and wait for a keystroke
; OUTPUT
;   outfile - name of FITS file name for joined spectrum
;             NOTE - if outfile is not specified, file2 assumed to be outfile

; set and check keywords

ocol=3
if (keyword_set(outfile) eq 0) then begin
  twoflag=0
  outfile=file2 
endif else twoflag=1     ; twoflag indicates that two input files are specified
if (keyword_set(bump) eq 0) then bumpflag=0 else bumpflag=1
if (keyword_set(bump2) eq 0) then bump2=0
if (keyword_set(nobump) eq 0) then nobumpflag=0 else nobumpflag=1
if (keyword_set(cross) eq 0) then crossflag=0 else crossflag=1

; read spectral FITS files, get number of segments, min, max order in each
; check to see that both are valid
; if sp1 is invalid, need to write sp2 data to sp1 and reset twoflag=0

sp1=readfits(file1,hdr1,/silent)
if (sp1[0] ne -1) then begin ; continue only if the file is valid
  nseg1=sxpar(hdr1,'NSEG')
  min1=min(sp1[ocol,*]) & max1=max(sp1[ocol,*])
endif

; read second file if there is one set (valid or not)

if (twoflag eq 1) then begin
  sp2=readfits(file2,hdr2,/silent)
  if (sp2[0] ne -1) then begin ; continue only if the file is valid
    nseg2=sxpar(hdr2,'NSEG')
    min2=min(sp2[ocol,*]) & max2=max(sp2[ocol,*])
  endif 
endif else sp2=-1           ; set sp2 for invalid file

; now branch for following conditions
;   both sp1 and sp2 are valid - checking to bump orders in sp2 is automatic
;   only sp1 is valid - bump only if bump keyword set
;   only sp2 is valid - 
;   both sp1 and sp2 are invalid - this block excludes plotting and writing
; if two input files, check for overlapping order numbers
; if necessary, increment order numbers in sp2
; if one input file and bump is set, increment order numbers in sp1

if (sp1[0] ne -1 and sp2[0] ne -1) then begin    ; both valid

; if orders overlap and nobumpflag not set, del = biggest possible delta

  if (max1 ge min2 and min1 le max2 and nobumpflag eq 0) then $
    del = max( [max1-min2+1,nseg1,nseg2] ) else del=0  

; if bump2 is set, then manually fix del = bump2 (regardless of order overlap)

   if (bump2 ne 0) then del = bump2

;  bump order number (or not, if del=0)

   sp2[ocol,*]=sp2[ocol,*]+del

;  check number of columns - add more columns to whichever array is smaller

   ncol1 = n_elements(sp1[*,0])
   ncol2 = n_elements(sp2[*,0])
   if (ncol1 ne ncol2) then begin
     if (ncol1 gt ncol2) then begin
       new = fltarr(ncol1,n_elements(sp2[0,*]))
       new[0:ncol2-1,*] = sp2
       sp2 = new
     endif else begin
       new = fltarr(ncol2,n_elements(sp1[0,*]))
       new[0:ncol1-1,*] = sp1
       sp1 = new
     endelse
   endif

;  concatenate spectral arrays

   sp=[[sp1],[sp2]]

;  use hdr1, but if both files valid and crossflag set, call crossradec

   if (crossflag eq 1) then outhdr=crossradec(hdr1,hdr2) else outhdr=hdr1

endif

if (sp1[0] ne -1 and sp2[0] eq -1) then begin    ; only sp1 is valid
  if (twoflag eq 0 and bumpflag eq 1) then begin ; bump if bump set
    del = max( [max1-min1+1,max1,nseg1] )        ; was max1-min1+1
    sp1[ocol,*]=sp1[ocol,*]+del                  ; bump order number
  endif
  sp=sp1
  outhdr=hdr1
endif

if (sp1[0] eq -1 and sp2[0] ne -1) then begin    ; only sp2 is valid

  if (nobumpflag eq 0) then begin
    del = max( [max2-min2+1,nseg2] )             ; was max2-min2+1
    sp2[ocol,*]=sp2[ocol,*]+del                  ; bump order number
  endif
  sp=sp2
  outhdr=hdr2
endif

if (sp1[0] eq -1 and sp2[0] eq -1) then begin    ; both files invalid
  print,'Warning - no valid data, cannot write file ',outfile ; warn and exit
endif else begin

; plot if requested

  if (keyword_set(plot) ne 0) then begin
    spplot,sp,title=outfile,_extra=e
    speplot,sp,_extra=e
    zq=get_kbrd(1)
  endif

; write output file
  
  wr_spfits,outfile,sp,-1,fits_hdr=outhdr

endelse

END
