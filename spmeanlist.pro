FUNCTION spmeanlist,splist,nspec,NORM=norm,MEDIAN=median

; 12 Jul 10 created
;
; spmeanlist finds the mean or median at each wavelength for a list of spectra
; and returns the resulting spectrum
; error column is the larger of:
;   1) smaller of standard deviation of normalized and unnormalized spectra
;      (divided by sqrt(nspec))
;   2) propagated errors
;
; INPUT
;   splist   - list of spectral data arrays to be combined
;              these are passed concatenated and thus must have the same length
;              it is assumed they have the same wavelength grid
;   nspec    - number of spectra in list (since no. of columns can vary)
;              if nspec unspeficied, assumed that spectra have 4 cols
;   median   - if set, returns median instead of mean
;   norm     - optional keyword to normalize spectra before finding median
;              DEFAULT is to simply find mean or median with no normalization
;              (normalization has no affect on mean)
; OUTPUT     - returns mean (or median) spectrum

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check array dimensions, load data into newsp, flux and error arrays

if (keyword_set(nspec) eq 0) then begin
  ncol=4
  nspec=n_elements(splist[*,0])/ncol
endif else begin
  ncol=n_elements(splist[*,0])/nspec 
endelse
nrow=n_elements(splist[0,*])

newsp         = fltarr(ncol,nrow)
newsp[lcol,*] = splist[lcol,*]
newsp[ocol,*] = splist[ocol,*]

flux   = fltarr(nspec,nrow)
error  = fltarr(nspec,nrow)
for i=0,nspec-1 do flux[i,*]  = splist[fcol+i*ncol,*]
for i=0,nspec-1 do error[i,*] = splist[ecol+i*ncol,*]

; normalize input spectra, flux and errors if norm keyword set

if (keyword_set(norm) ne 0) then begin
  meanval=fltarr(nspec)
  for i=0,nspec-1 do meanval[i] = mean(flux[i,*],/nan) ; mean of each spectrum
  meanval=meanval/mean(meanval)                        ; normalize meanval to 1
  for i=0,nspec-1 do flux[i,*]=flux[i,*]/meanval[i]    ; normalize spectra
  for i=0,nspec-1 do error[i,*]=error[i,*]/meanval[i]  ; normalize errors
  print,meanval
endif

; compute new fluxes

if (keyword_set(median) eq 0) then medflag=0 else medflag=1
for j=0,nrow-1 do begin
  if (medflag eq 0) then newsp[fcol,j]=mean(flux[*,j],/nan) $
    else newsp[fcol,j]=median(flux[*,j],/even)
endfor

; find error - the complicated part

err_nonorm = fltarr(nrow) ; uncertainty of unnormalized spectra
err_norm   = fltarr(nrow) ; uncertainty of normalized spectra
err_prop   = fltarr(nrow) ; propagated error

for j=0,nrow-1 do begin

; error without normalizing input spectra

  err_nonorm[j] = stddev(flux[*,j])/sqrt(nspec)

; propagated error - two factors of sqrt(nspec), once to avg rms, once for unc.

  err_prop[j]   = sqrt(total(error[*,j]^2))/nspec 

endfor

; if input spectra not already normalized, then do it now

if (keyword_set(norm) eq 0) then begin
  meanval=fltarr(nspec)
  for i=0,nspec-1 do meanval[i] = mean(flux[i,*],/nan) ; mean of each spectrum
  meanval=meanval/mean(meanval)                        ; normalize meanval to 1
  for i=0,nspec-1 do flux[i,*]=flux[i,*]/meanval[i]    ; normalize spectra
endif

; error after normalizing input spectra

for j=0,nrow-1 do begin

  err_norm[j] = stddev(flux[*,j])/sqrt(nspec)

; select smaller of err_norm and err_nonorm
; and larger of that result and err_prop

  value=min([err_norm[j],err_nonorm[j]],/nan)
  newsp[ecol,j] = max([value,err_prop[j]],/nan)

endfor

RETURN,newsp
END
