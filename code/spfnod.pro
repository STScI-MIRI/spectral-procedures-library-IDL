PRO spfnod,file1,file2,outfile,NOPAIR=nopair,PLOT=plot,_extra=e

; 24 Apr 07 created by modifying spfnodjoin to process only 2 spectra, not 4
;
; spfnod combines two spectra on the same wavelength grid
; they must have the same number of columns
; the header from the first file will be used as the output FITS header
; assumed that col 3 = order number
; nods are combined using sppair.pro (or spmean.pro if NOPAIR set)
;
; INPUT
;   file1,2 - spectral FITS files for the two spectra to combine 
;   outfile - name of output spectral FITS file
;   nopair  - optional keyword to simply average spectra for each segment
; OUTPUT - writes outfile to disk

ocol=3

sp1=readfits(file1,hdr,/silent)
sp2=readfits(file2,/silent)

if (n_elements(nopair) eq 0) then begin
  spectrum=sppair(sp1,sp2,_extra=e)
endif else begin
  spectrum=spmean(sp1,sp2,_extra=e)
endelse

m0=min(spectrum[ocol,*]) & m1=max(spectrum[ocol,*])
for m=m0,m1 do begin
  o_idx=where(spectrum[ocol,*] eq m)
  if (max(o_idx) gt -1) then begin
    if (m eq m0) then olen=n_elements(o_idx) $
    else              olen=[olen,n_elements(o_idx)]
  endif
endfor

if (keyword_set(plot) ne 0) then begin
  spplot,spectrum,title=outfile,_extra=e
  speplot,spectrum
  zq=get_kbrd(1)
endif

wr_spfits,outfile,spectrum,olen,fits_hdr=hdr

END
