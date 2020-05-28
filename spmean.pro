FUNCTION spmean,sp1,sp2,NONORM=nonorm,NOCHECK=nocheck,GRIDFILE=gridfile

;  7 Dec 11 added gridfile keyword to pass wavelength grid for regridding
; 26 Feb 08 updated error calculation if S/N > 1000
; 12 Nov 05 returns error in first spectrum if error=0
; 15 Mar 05 added the nonorm and nocheck keywords 
;           the default is now to find errors before and after normalizing nods
;             and to use the smaller mean errors
;           nocheck means to use the normalization method
;           nonorm means to use the non-normalized nods to find the error
; 17 Dec 04 created
;
; spmean finds the average at each wavelength between sp1 and sp2
; and returns the resulting spectrum
; error column is the uncertainty in the mean (which is half the diff)
;   errors calculated after segments are normalized to eliminate photometric
;   effects
;
; INPUT
;   sp1,sp2  - spectral data arrays to be combined
;              if they have different lengths, sp1 is returned
;              it is assumed they have the same grid
;   gridfile - keyword to regrid; argument is the wavelength FITS file 
;   nonorm   - keyword to stop normalization of nods when finding error
;   nocheck  - keyword to use normalized nods to find errors
;   (default)  use whichever method finds smaller errors
;              NOTE - if both nonorm and nocheck are set, nonorm is used
; OUTPUT - returns mean spectrum

lcol=0 & fcol=1 & ecol=2 & ocol=3

; if gridfile keyword set, load wavelength grid and regrid

if (keyword_set(gridfile) ne 0) then begin
  lamgrid=readfits(gridfile,/silent)
  sp1=spregrid(sp1,lamgrid)
  sp2=spregrid(sp1,lamgrid)
endif

len1=n_elements(sp1[lcol,*])
len2=n_elements(sp2[lcol,*])

if (len1 eq len2) then begin

  newsp=sp1
  newsp[fcol,*]=0.5*(sp1[fcol,*]+sp2[fcol,*])

; find error - the complicated part

;   error without normalizing nods

    err_nonorm=0.5*abs(sp1[fcol,*] - sp2[fcol,*])

;   error with normalization of nods

    fm1=median(sp1[fcol,*])
    fm2=median(sp2[fcol,*])
    fmavg=0.5*(fm1+fm2)     ; used to normalize fluxes before calculating errors
    err_norm=0.5*abs(sp1[fcol,*]*fmavg/fm1 - sp2[fcol,*]*fmavg/fm2)

    if (n_elements(nonorm) eq 1) then begin
      newsp[ecol,*]=err_nonorm
    endif else begin
      if (n_elements(nocheck) eq 1) then begin
        newsp[ecol,*]=err_norm
      endif else begin
        if (mean(err_nonorm) lt mean(err_norm)) then newsp[ecol,*]=err_nonorm $
        else newsp[ecol,*]=err_norm
      endelse
    endelse

; if error < 0.1% flux, then replace with max of it or error from first spec.

    idx=where(newsp[ecol,*] lt 0.001*abs(newsp[fcol,*]),icnt)
    if (icnt gt 0) then for i=0,icnt-1 do begin
      newsp[ecol,idx[i]]=max([sp1[ecol,idx[i]],newsp[ecol,idx[i]]])
    endfor

endif else begin

  print,'Warning in spmean.  Spectra have different lengths, returning first.'
  newsp=sp1

endelse

RETURN,newsp
END
