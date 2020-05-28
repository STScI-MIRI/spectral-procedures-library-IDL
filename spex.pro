FUNCTION spex,spin,spcin,l0,l1,lam_center,sigratio,ORDER=order,MODE=mode,FLUX=flux,RATIO=ratio,DIAG=diag,PLOT=plot,FLAM=flam,_extra=e

; 11 Mar 15 fixed a resizing bug when passed an array that doesn't have 4 cols
; 14 Feb 13 added flam keyword to assume input is in F_lam units (um^-1)
;           which effectively turns off the conversion to F_nu
; 16 Sep 10 modified to handle arrays smaller than needed
;  3 Jul 08 added _extra keyword (belatedly!)
;  2 Jul 08 added plotting keyword
; 27 Feb 07 moved all error propagation from spcex to spex
;           must now have sigratio - uncertainty in continuum fit
; 14 Mar 06 shifted lam_center to red by one-half wavelength interval
;           tested this with various input spectra, now getting correct center
; 12 Sep 05 added mode keyword as an alternative way of invoking flux or ratio 
;           also allowing spcin=scalar zero now (for no background)
;  6 Sep 05 only sorting if order not set, sorting dellam vector properly now
; 31 Aug 05 renamed as spex, added error to lam_center, /ratio keyword
; 30 Aug 05 added lam_center parameter to find central wavelength
; 29 Mar 05 added /flux keyword
;  5 Mar 05 modified to use internal uncertainty when available
; 15 Oct 04 created as spexw
;
; given a spectral data array and corresponding continuum
; determine equivalent width of feature from l0 to l1 (or flux or flux ratio)
; both spectra data arrays must be on the same wavelength grid
;
; assumed that flux is in Jy
; flux and error units converted to W m^-2 um^-1
;
; INPUT
;   spin       - spectral data array containing the original spectrum
;   spcin      - spectral data array containing the continuum estimate
;              - note, if spcin = scalar 0, then a vector of zeroes is used
;   l0,l1      - starting and ending wavelength to find equivalent width
;   order      - keyword to limit analysis to one order
;   mode       - optional parameter to set flux (1) or ratio (2) keywords
;   flux       - keyword to return equiv. flux instead of equiv. width (mode=1)
;   ratio      - keyword to return ratio of eq flux to continuum (mode=2)
;   lam_center - parameter used to return central wavelength and error
;                defined as wavelength with half the flux on either side
;   sigratio   - optional parameter to pass fractional error in continuum flux
;   diag       - optional keyword to turn on diagnostic mode
;   plot       - optional keyword to plot the fit
;
; OUTPUT - returns equivalent width and uncertainty as a two-element vector

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check dimensions of spin, expand if necessary, copy into spdata

if (n_elements(spin[*,0]) lt 4) then begin
  spdata=fltarr(4,n_elements(spin[0,*]))
  spdata[0:n_elements(spin[*,0])-1,*] = spin
endif else spdata=spin

; check keywords - mode supersedes flux or ratio keywords

if (keyword_set(mode) eq 0) then begin
  mode=0
  if (keyword_set(flux) eq 1) then mode=1
  if (keyword_set(ratio) eq 1) then mode=2
endif 

if (n_elements(order) gt 0) then ordflag=1 else ordflag=0

; if spcin = scalar 0, then load spcont with vector of zeroes
; if spcin has too few columns, expand as spcont loaded

if (n_elements(spcin) eq 1 and total(spcin) eq 0) then begin
  spcont=spdata & spcont[fcol,*]=0 & spcont[ecol,*]=0
endif else if (n_elements(spcin[*,0]) lt 4) then begin
  spcont=fltarr(4,n_elements(spcin[0,*]))
  spcont[0:n_elements(spcin[*,0])-1,*] = spcin
endif else spcont=spcin

; check size of arrays

if (n_elements(spdata[lcol,*]) ne n_elements(spcont[lcol,*])) then begin
  print,'Error in spex.  Spectrum and continuum have different lengths.'
endif

; check order keyword and either load relevant data into spo and spc or sort

spo=spdata & spc=spcont 
if (ordflag eq 1) then begin
  idx=where(spdata[ocol,*] eq order)
  if (max(idx) gt -1) then begin
    spo=spdata[*,idx]
    spc=spcont[*,idx]
  endif else begin
    print,'Warning in spex.  No data in order, using full spectrum.'
    ordflag=0
  endelse
endif 

if (ordflag eq 0) then begin
  idx=sort(spo[lcol,*])
  spo=spo[*,idx]
  spc=spc[*,idx]
endif

; set up dellam vector - contains delta lambda for each pixel

dellam=abs(reform(shift(spo[lcol,*],1)-spo[lcol,*]))
dellam[0]=dellam[1]                        ; must replace first pixel 

; convert spectral data from F_nu units to F_lambda units 

if (keyword_set(flam) eq 0) then begin
  spo[fcol,*] = spo[fcol,*] * 3e-12 / (spo[lcol,*]^2)
  spo[ecol,*] = spo[ecol,*] * 3e-12 / (spo[lcol,*]^2)
  spc[fcol,*] = spc[fcol,*] * 3e-12 / (spc[lcol,*]^2)
endif

; find wavelengths for summing equivalent width

idx=where(spo[lcol,*] ge l0 and spo[lcol,*] le l1)

if (max(idx) gt -1) then begin

; NOW finding S/N ratio from error column
; had been using smoothing flux in this range 
; also must find mean wavelength element
; use these to find uncertainty in equivalent width

;  spsmoo=spsmooth(spsmooth(spo[*,idx],7),7)
;  snr=mean(spo[fcol,idx])/stddev(spsmoo[fcol,*]-spo[fcol,idx])
  snr=mean(spo[fcol,idx]/spo[ecol,idx])
  wavel=mean(dellam)
  sigw=sqrt(2.0*n_elements(idx))*wavel/snr
  sigf=sqrt(total(spo[ecol,idx]^2*dellam[idx]^2))

; find equivalent width, eq flux, integrated flux of continuum

  eqw=total( (1.0-(spo[fcol,idx]/spc[fcol,idx])) * dellam[idx] )
  eqf=total((spo[fcol,idx]-spc[fcol,idx])*dellam[idx])
  eqc=total(spc[fcol,idx]*dellam[idx])

; return value depends on mode, 0 - eq width, 1 - eq flux, 2 - flux ratio

  case mode of
    0 : retval=[eqw,sigw]
    1 : retval=[eqf,sigf]
    2 : retval=[eqf/eqc,sigf/eqc] ; divide flux, error by integrated continuum
  endcase

; propagate uncertainty in continuum into returned uncertainty

; convert sigratio, then add in quadrature to original uncertainty
; if finding eq width, convert to um by multiplying by delta lambda (l2-l1)
; if finding eq flux, multiply by mean integrated continuum (W/m^2)
; if finding flux ratio, use unmodified sigratio
; after all that, then find new rms error

  if (keyword_set(sigratio) eq 0) then sigratio=0.0 ; set if not supplied
  case mode of
    0 : unc_cont=sigratio*(l1-l0)
    1 : unc_cont=sigratio*eqc
    2 : unc_cont=sigratio
  endcase

  retval[1]=sqrt(retval[1]^2+unc_cont^2)

; find lam_center - central wavelength of feature
; to find error, find wavelength corresponding to half flux +/- sigval and avg
; sigval is defined to be the fractional uncertainty in the feature strength
;   with 27 Feb 07 change, this now includes propagated uncertainty in continuum

  sumf=(spo[fcol,idx]-spc[fcol,idx])*dellam[idx]
  for j=n_elements(idx)-1,0,-1 do sumf[j]=total(sumf[0:j])
  sumf=sumf/eqf

; first term below is just fractional random error in extracted feature flux
; second term is the error in continuum flux / feature flux

  sigval=sqrt((sigf/eqf)^2 + (sigratio*eqc/eqf)^2)

; now that we have sigval and summed flux, loope
; loop 0 - fit to half-sigval, 1 - to half, 2 - to half+sigval
; lam_center= result from loop 1, error = mean diff from loop1 to loops 0 and 2

  lam_c=fltarr(3)
  spot=[0.5-sigval,0.5,0.5+sigval]
  for i=0,2 do begin
    i0=max(where(sumf le spot[i])) & if (i0 le -1) then i0=0
    i1=min(where(sumf ge spot[i])) & if (i1 le -1) then i1=n_elements(idx)-1
    lam0=spo[lcol,idx[i0]]+0.5*dellam[idx[i0]]
    lam1=spo[lcol,idx[i1]]+0.5*dellam[idx[i1]]
    if (lam0 eq lam1) then lam_c[i]=lam0 else begin
      lam_c[i] = lam0 + (lam1-lam0)*(spot[i]-sumf[i0]) / (sumf[i1]-sumf[i0])
      if (keyword_set(diag) eq 1) then print,i0,lam0,sumf[i0]
      if (keyword_set(diag) eq 1) then print,i1,lam1,sumf[i1]
    endelse
  endfor
  lam_center=[lam_c[1],0.5*abs(lam_c[2]-lam_c[0])]

endif else begin

  print,'Warning in spex.  No data in range, returning zeroes.'
  retval=[0.0,0.0]
  lam_center=[0.0,0.0]

endelse

if (keyword_set(plot) eq 1) then begin
  spplot,spdata,_extra=e
  spplot,spcont,_extra=e,/over
  ymin=spmin(spdata) & ymax=spmax(spdata)
  if (ymin gt 0) then ymin=0.5*ymin else ymin=2.0*ymin
  if (ymax gt 0) then ymax=2.0*ymax else ymax=0.5*ymax
  oplot,[l0,l0],[ymin,ymax],_extra=e
  oplot,[l1,l1],[ymin,ymax],_extra=e
  zq=get_kbrd(1)
endif

return,retval
END
