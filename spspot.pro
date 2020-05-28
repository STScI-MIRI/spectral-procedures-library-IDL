FUNCTION spspot,sp,lam,RJ=rj,ORDER=order,ERROR=error,SILENT=silent

; 28 Jun 13 added silent keyword
; 25 Jun 13 added error keyword to extract from the error col instead of flux
;  3 Jan 11 fixed a bug - order specification was ignored if wavelength matched
; 29 Apr 10 improved how to handle case of no data in range
; 21 Jun 06 added /rj keyword
;  4 Nov 04 ensure that returned value is not an array of one element
; 14 Oct 04 created
;
; spspot determines the flux in a spectrum at a given wavelength
; if the wavelength appears in the spectrum exactly, the flux is returned
; else the interpolated flux between the two neighbors is returned
; if lam out of range of spectrum, min or max is returned with a warning
;
; INPUT
;   sp     - spectral data array
;   lam    - wavelength in units matching those in col 0 of spectral data array
;   order  - keyword to limit analysis to one order
;   rj     - keyword to return value in Rayleigh-Jeans units
;   silent - keyword to repress warnings
; OUTPUT - returns flux at lam, exact or interpolated

lcol=0 & fcol=1 & ocol=3
spo=sp

if (keyword_set(error) ne 0) then fcol=2 ; point to error col instead of flux

if (n_elements(order) eq 1) then begin
  idx=where(sp[ocol,*] eq order)
  if (max(idx) gt -1) then spo=sp[*,idx] else begin
    if (keyword_set(silent) ne 0) then $
      print,'Warning in spspot.  No data in order, using full array'
   endelse
endif

idx=where(spo[lcol,*] eq lam)

if (max(idx) gt -1) then begin ; matched wavelength exactly

  retval=mean(spo[fcol,idx])

endif else begin ; either out of range or must interpolate

  lessidx=where(spo[lcol,*] lt lam)
  moreidx=where(spo[lcol,*] gt lam)

  if (max(lessidx) gt -1) then begin ; have data to blue
    lessflag=1
    spless=spo[*,lessidx]
    maxless=max(spless[lcol,*],spot0)
  endif else lessflag=0
  if (max(moreidx) gt -1) then begin ; have data to red
    moreflag=1
    spmore=spo[*,moreidx]
    minmore=min(spmore[lcol,*],spot1)
  endif else moreflag=0

  if (lessflag eq 1 and moreflag eq 1) then begin ; interpolate
    l0=spless[lcol,spot0] & l1=spmore[lcol,spot1]
    f0=spless[fcol,spot0] & f1=spmore[fcol,spot1]
    retval=f0 + (f1-f0)*(lam-l0)/(l1-l0)
  endif else begin
    if (keyword_set(silent) eq 0) then $
      print,'Warning in spspot.  Data out of range, returning closest value.'
    if (lessflag eq 1) then begin ; no data to red of lam, take max value
      maxdata=max(spo[fcol,*],maxspot)
      retval=maxdata
    endif
    if (moreflag eq 1) then begin ; no data to blue of lam, take min value
      mindata=min(spo[fcol,*],minspot)
      retval=mindata
    endif
  endelse  

endelse

if (keyword_set(rj) eq 0) then return,retval else return,retval*lam^2
END
