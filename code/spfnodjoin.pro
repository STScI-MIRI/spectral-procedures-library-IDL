PRO spfnodjoin,file1a,file1b,file2a,file2b,outfile,NOPAIR=nopair,PLOT=plot,_extra=e

; 20 Aug 05 added _extra to spplot call
; 15 Mar 05 added _extra to spmean call
; 29 Dec 04 added a get_kbrd request after a plot
; 18 Dec 04 added PLOT keyword for plotting; default is now not to plot
; 17 Dec 04 added NOPAIR keyword to bypass sppair and simply average nods
;  8 Mar 04 copied and modified to make call to sppair for nod pos'ns
;  8 Mar 04 renamed from spjoin to spfjoin
;  2 Feb 04 pass header to wr_spfits, overlooked this before
; 27 Dec 03 created as spjoin
;
; spfnodjoin combines four spectral FITS files (two nods, two orders)
; they must have the same number of columns
; the header from the first file will be used as the output FITS header
; if order numbers overlap, then order numbers from file 2 are incremented
; assumed that col 3 = order number
; nods are combined using sppair.pro (or spmean.pro if NOPAIR set)
;
; INPUT
;   file1a,1b - spectral FITS files for two spectra to combine in first segment
;   file2a,2b - ditto for second segment
;   outfile   - name of output spectral FITS file
;   nopair    - optional keyword to simply average spectra for each segment
; OUTPUT - writes outfile to disk

sp1a=readfits(file1a,hdr,/silent)
sp1b=readfits(file1b,/silent)
sp2a=readfits(file2a,/silent)
sp2b=readfits(file2b,/silent)

if (n_elements(nopair) eq 0) then begin
  sp1=sppair(sp1a,sp1b,_extra=e)
  sp2=sppair(sp2a,sp2b,_extra=e)
endif else begin
  sp1=spmean(sp1a,sp1b,_extra=e)
  sp2=spmean(sp2a,sp2b,_extra=e)
endelse

; from this point on, the procedure is the same as spfjoin.pro

; check for overlapping order numbers
; if necessary, increment order numbers in sp2

min1=min(sp1[3,*]) & max1=max(sp1[3,*])
min2=min(sp2[3,*]) & max2=max(sp1[3,*])

overflag=0
if ((max1 ge min2) and (min1 le max2)) then overflag=1

if (overflag eq 1) then begin
  del=max1-min2+1
  sp2[3,*]=sp2[3,*]+del
endif

; concatenate spectral arrays

sp=[[sp1],[sp2]]

; determine order length array needed by wr_spfits

m0=min(sp[3,*]) & m1=max(sp[3,*])
for m=m0,m1 do begin
  o_idx=where(sp[3,*] eq m)
  if (max(o_idx) gt -1) then begin
    if (m eq m0) then olen=n_elements(o_idx) $
    else              olen=[olen,n_elements(o_idx)]
  endif
endfor

if (keyword_set(plot) ne 0) then begin
  spplot,sp,title=outfile,_extra=e
  speplot,sp
  zq=get_kbrd(1)
endif

wr_spfits,outfile,sp,olen,fits_hdr=hdr

END
