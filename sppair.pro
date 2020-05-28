FUNCTION sppair,sp1,sp2,SIGMA=sigma,WIDTH=width,GRIDFILE=gridfile,SPECERROR=specerror,DIAG=diag,BOB=bob

;  2 Mar 12 added bob keyword to remove one pixel from each end of each order
; 13 Dec 11 clean-up of code and comments - no longer modifying input spectra
;  7 Dec 11 added the gridfile keyword to pass wavelength grid for regridding
;  3 Jun 10 added the specerror option
; 26 Feb 08 updated error calculation if S/N > 1000
; 17 Mar 07 preventing oscillations between spectra when they differ
; 12 Nov 05 if error near zero, return error in first spectrum
; 19 Dec 04 repaired bug with weighting
; 15 Dec 04 reduce std dev to uncertainty in mean where both nods good
; 21 Mar 04 either reject or use data, do not use continuum, return uncertainty
;  8 Mar 04 created
;
; given two spectral data arrays, combine them and eliminate spikes
; algorithm:
;   generate continuum spectra for each input spectrum
;   generate difference spectra for each = abs(continuum-spectrum)
;   generate pass spectrum:
;   (1) at each wavelength, choose min of the two diff spectra
;   (2) replace small pass values by raising to median or mean of pass array
;   generate a weight array by comparing differences to pass array
;     w=0 if diff > pass*sigma, otherwise w=1
;     if both weights=0 for a given wavelength element, reset to 1
;     but don't set w=0 if the spectra vary more than the spike varies
;   combine spectra by averaging
;   compute uncertainty as though all data were used
; INPUT
;   sp1,sp2   - spectral data arrays, copied to spa,spb, which are modified
;   sigma     - number of times pass before data replaced
;   sigma     - keyword, multiple of pass array to exceed for data replacement
;   width     - keyword - size of box for median calls
;   gridfile  - keyword to regrid; argument is the wavelength FITS file
;   specerror - keyword to calculate errors after normalizing spectra
;             - this removes the systematic (photometric) error
;   diag      - keyword to turn diagnostic plotting on
; OUTPUT
;   spike-corrected and averaged spectrum

; setup

diagmode=2
lcol=0 & fcol=1 & ecol=2 & ocol=3

; check keywords

if (keyword_set(sigma) eq 0) then sigma=5.0
if (keyword_set(width) eq 0) then width=5
if (keyword_set(diag) eq 0)  then diagflag=0 else diagflag=diag

; if bob keyword set, then call spbob to remove a pixel from ends of each order

if (keyword_set(bob) ne 0) then begin
  spa=spbob(sp1)
  spb=spbob(sp2)
endif else begin
  spa=sp1
  spb=sp2
endelse

; if gridfile keyword set, load wavelength grid and regrid

if (keyword_set(gridfile) ne 0) then begin
  lamgrid=readfits(gridfile,/silent)
  spdum=sporegrid(spa,lamgrid)
  spa=spdum
  spdum=sporegrid(spb,lamgrid)
  spb=spdum
endif 

; check that spectral arrays are the same size
; and initialize new spectral data array newspec

len=n_elements(spa[lcol,*])
if (len ne n_elements(spb[lcol,*])) then begin
  print,"Error.  Spectral sizes do not match."
  stop
endif
newspec=spa

; load wavelength, flux, continuum, diff spectra, pass spectrum
; continuum spectra are media of flux spectra

lam   = reform(spa[lcol,*])
fluxa = reform(spa[fcol,*])
fluxb = reform(spb[fcol,*])
erra  = reform(spa[ecol,*])
cona  = median(fluxa,width)
conb  = median(fluxb,width)
diffa = abs(fluxa-cona)
diffb = abs(fluxb-conb)
diffx = 0.5*abs(cona-conb) 
pass  = fltarr(len)

; generate pass array, each pixel is min of diffa or diffb

for i=0,len-1 do pass[i]=min([diffa[i],diffb[i]]) ; minimum difference

if (diagflag eq 1) then begin
  spplot,spa
  spplot,spb,/over
endif
if (diagflag eq 2) then begin
   plot,lam,diffa
   oplot,lam,diffb
endif

; raise low values of pass to median of pass array 

;minpass=max([median(pass),mean(pass)]) ; used to be smaller of median or mean

minpass=median(pass)
bad_idx=where(pass lt minpass)
if (diagflag gt 0) then $
  print,'median ',median(pass),' mean ',mean(pass),' minpass ',minpass
if (bad_idx[0] gt -1) then pass[bad_idx] = minpass
if (diagflag eq 2) then oplot,lam,sigma*pass,th=2

; set weight arrays wt_a and wt_b
; default is 1.0
; set weights to zero when two conditions fulfilled
;   1) pixel is a spike (diff > sigma*pass)
;   2) pixel deviates more than the one in the other nod
; if both wt_a and wt_b for a given wavelength element=0, reset to 1
;   this shouldn't happen, but it would indicate the pixel is hopeless

wt_a=fltarr(len)+1.0
wt_b=wt_a

bad_idx=where(diffa gt sigma*pass and diffa gt diffx)
if (bad_idx[0] gt -1) then begin
  wt_a[bad_idx]=0.0
  if (diagflag eq 1) then oplot,lam[bad_idx],fluxa[bad_idx],psym=4
  if (diagflag eq 2) then oplot,lam[bad_idx],diffa[bad_idx],psym=4
endif

bad_idx=where(diffb gt sigma*pass and diffb gt diffx)
if (bad_idx[0] gt -1) then begin
  wt_b[bad_idx]=0.0
  if (diagflag eq 1) then oplot,lam[bad_idx],fluxb[bad_idx],psym=6
  if (diagflag eq 2) then oplot,lam[bad_idx],diffb[bad_idx],psym=6
endif

bad_idx=where(wt_a eq 0.0 and wt_b eq 0.0)
if (bad_idx[0] gt -1) then begin
  print,'There are ',n_elements(bad_idx),' hopeless pixels'
  wt_a[bad_idx] =1.0
  wt_b[bad_idx] =1.0
endif

; this block added 3 Jun 10
; if specerror set, normalize fluxb to fluxa - save as normb
; use normb in place of fluxb to determine error in the next step

normb=fluxb ; set normb=fluxb first - default if specerror not set
if (keyword_set(specerror) ne 0) then begin
  idx=where(finite(fluxa) and finite(fluxb))
  if (idx[0] ne -1) then normb=normb*(mean(fluxa[idx])/mean(fluxb[idx]))
endif 

; combine spectra by weighted averaging
; generate new uncertainties
; divide uncertainties by sqrt(2) where both nods good

newspec[fcol,*]=(fluxa*wt_a+fluxb*wt_b)/(wt_a+wt_b)
for i=0,len-1 do newspec[ecol,i]=stddev([fluxa[i],normb[i]])
idx=where(newspec[fcol,*] ne fluxa and newspec[fcol,*] ne fluxb)
if (max(idx) gt -1) then newspec[ecol,idx]=newspec[ecol,idx]/sqrt(2.0)

; if error < 0.1% flux, then replace with max of it or error from first spec.

idx=where(newspec[ecol,*] lt 0.001*abs(newspec[fcol,*]),icnt)
if (icnt gt 0) then for i=0,icnt-1 do begin
  newspec[ecol,idx[i]]=max([erra[idx[i]],newspec[ecol,idx[i]]])
endfor

if (diagflag eq 1) then spplot,newspec,th=2,/over

return,newspec

END
