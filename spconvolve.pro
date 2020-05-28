FUNCTION spconvolve,spin,fwhm,ORDER=order,RANGE=range

; 28 Jun 13 created
;
; given an input spectrum, convolve it with a Gaussian of FWHM fwhm
; and return the result
;
; NOTE - need to PROPAGATE the error better to reflect the smoothing
;
; INPUT
;   spin  - input spectral data array
;   fwhm  - FWHM in um (not pixels) (= 2.355 sigma) - NO DEFAULT - must be given
;   order - keyword to limit the convolution to an order
;   range - keyword to limit integration to N sigma, default=3
; OUTPUT

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check keywords

if (keyword_set(fwhm) eq 0) then begin
  print,'Error.  Must set FWHM.  Stopping.'
  stop
endif

; determine sigma from FWHM and redefine range in um

sigma=fwhm/2.355
if (keyword_set(range) eq 0) then range=3*sigma else range=range*sigma

; copy needed data to sp and if needed, limit it to the order requested

goflag=1
sp=spin

if (keyword_set(order) ne 0) then begin
  idx=where(sp[ocol,*] eq order)
  if (idx[0] gt -1) then sp = spin[*,idx] else begin
    print,'Warning.  No data in requested order, returning input spectrum.'
    goflag=0
  endelse
endif 

spout=sp

if (goflag eq 1) then begin

  len=n_elements(spout[lcol,*])

  for i=0,len-1 do begin

; find all elements within range of gridpoint

    lam0=spout[lcol,i]
    idx = where(sp[lcol,*] ge lam0-range and sp[lcol,*] le lam0+range)

; nothing in range (shouldn't happen), just use input wavelength element

    if (idx[0] lt 0) then begin

      spout[fcol,i] = sp[fcol,i] 
      spout[ecol,i] = sp[ecol,i] 

    endif else begin

; compute gaussian, sum of weights, weighted flux, etc.

      wt    = exp(-(sp[lcol,idx]-lam0)^2/(2*sigma^2))
      wtsum = total(wt)
      fsum = total(wt*sp[fcol,idx])
      esum = total(wt*sp[ecol,idx]) ; should RENORMALIZE to address smoothing

; find new flux and error - used input data if sum of wts = 0

      if (wtsum gt 0.0) then begin
        flux = fsum/wtsum  & error = esum/wtsum
      endif else begin
        flux = sp[fcol,i]  & error = sp[ecol,i]
      endelse

      spout[fcol,i] = flux & spout[ecol,i] = error

    endelse
  endfor

endif

RETURN,spout
END
