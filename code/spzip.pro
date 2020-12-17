FUNCTION spzip,insp1,insp2,SPECERROR=specerror,MODE=mode,_extra=e

; 15 Dec 11 may be done with the modifications
; 12 Dec 11 created by modifying code written by Vianney Lebouteiller
;
; original code:  /home/scripts/CASSIS_v4/cassis_combine_spectra.pro, line 335+
; or /home/scripts/CASSIS_v5/cassis_combine_spectra_fullsampling.pro
;
; combine two spectra on separate wavelength covering roughly the same range
; the code makes allowances for
; differences in slope or overall strength of the two spectra
;
; mmt_transform was written by J. Starck and is copied locally 
;   cf Starck et al. (1999, A&AS, 134, 135) 
; this code is considerably different from Vianney's code
; it uses the basic algorithm of a smoothed difference spectrum to adjust
;   the two input spectra before zipping them together (to avoid bad zig-zags)
; this code gives the option of using bspline_iterfit or mmt_transform to
;   smooth the difference spectrum
; and it uses a bad-pixel identification algorithm based on sppair.pro
;   and not the algorithm Vianney uses in his code
;
; OPEN ISSUES - should we regrid smoodiff before correcting sp2?
;
; INPUT
;   insp1,insp2 - the two spectra on different grids
;   mode        - determines how the difference spectrum will be smoothed
;               - 0 - DEFAULT - use bspline_iterfit, 1 - use mmt_transform
;   specerror   - keyword - to calculate errors after normalizing spectra
;               - this removes the systematic (photometric) error
; KEYWORDS PASSED TO SPSPIKE VIA _EXTRA IF DECLARED
;   sigma       - threshold for spike rejection (multiple of modified noise)
;   width       - smoothing box size for generation of median spectra
; OUTPUT        - the two spectra corrected and zipped into one sequence
;
; calls spspike.pro to remove spikes and divots

; setup - column definitions and dimensions, check keywords

lcol=0 & fcol=1 & ecol=2 & ocol=3 ; ocol not currently used

ncol = n_elements(insp1[*,0])
if (keyword_set(mode) eq 0) then mode=0

; set-up for iterating through orders and do it

minorder = min([min(insp1[ocol,*]),min(insp2[ocol,*])])
maxorder = max([max(insp1[ocol,*]),max(insp2[ocol,*])])

for m=minorder,maxorder do begin

  mdx1 = where(insp1[ocol,*] eq m) ; valid data in sp1
  mdx2 = where(insp2[ocol,*] eq m) ; valid data in sp2

  if (mdx1[0] gt -1 and mdx2[0] gt -1) then begin ; both spectra valid in range

;   set up arrays

    sp1 = insp1[*,mdx1]
    sp2 = insp2[*,mdx2]

    len1     = n_elements(sp1[0,*])
    len2     = n_elements(sp2[0,*])

    totalrow = len1+len2
    newerr1  = fltarr(len1)
    newerr2  = fltarr(len2)
    ordersp  = fltarr(ncol,totalrow)

;   copy spectral arrays into vectors

    lam1   = reform(sp1[lcol,*]) & lam2  = reform(sp2[lcol,*]) 
    flux1  = reform(sp1[fcol,*]) & flux2 = reform(sp2[fcol,*])
    err1   = reform(sp1[ecol,*]) & err2  = reform(sp2[ecol,*])

  ; determine the difference spectrum by regridding sp2 to sp1 and comparing

;    modsp2  = spregrid(sp2,sp1) ; spregrid does not appear to be robust enough
;    diff    = flux1 - modsp2[fcol,*]
    modflux2 = interpol(flux2,lam2,lam1)
    diff     = flux1 - modflux2

;   smooth the difference spectrum with either mmt_transform or bspline_iterfit

    if (mode eq 0) then begin ; use bspline_iterfit

;     determine break-point spacing for bspline_iterfit from the wavelength grid
;     spline-smooth the difference

      dlam   = median(lam1-shift(lam1,-1))
      breaks = 20*dlam                   ; 20 is somewhat arbitrary
      cc=bspline_iterfit(lam1,diff,bkspace=breaks,nord=3,yfit=smoodiff)

    endif else begin

;   nscales parameter depends on length of input to mmt_transform
;   then call mmt_transform
;   Vianney's trick is to pad the array with its reverse on both sides
;     to prevent ringing at the edges

      nscales = floor( alog10( 0.2*len1) / alog10(2.) )
      tt=mmt_transform([reverse(diff),diff,reverse(diff)],nscales)
      smoodiff=tt[nscales,len1:2*len1-1] ; pull last column out

    endelse

;   use smoodiff to correct the spectra
;   note that regridding for sp2 doesn't appear to be very important

    newflux1 = flux1 - 0.5*smoodiff
    newflux2 = flux2 + 0.5*interpol(smoodiff,lam1,lam2)

;   scale the errors multiplicatively by the shift in flux

    newerr1 = err1 * newflux2/flux1
    newerr2 = err2 * newflux2/flux2

;   generate new errors from the difference array and use where larger
;   if specerror is set, then use diff between new arrays (regridded)

    if (keyword_set(specerror) eq 0) then begin

      newerror1 = abs(diff)/2.0 
      newerror2 = newerror1

    endif else begin

      diff1  = newflux1 - interpol(newflux2,lam2,lam1)
      diff2  = newflux2 - interpol(newflux1,lam1,lam2)

      newerror1 = abs(diff1)/2.0
      newerror2 = abs(diff2)/2.0

    endelse

    idx=where(newerr1 lt newerror1)
    if (idx[0] gt -1) then newerr1[idx] = newerror1[idx]
    idx=where(newerr2 lt newerror2)
    if (idx[0] gt -1) then newerr2[idx] = newerror2[idx]

;   build new spectral arrays

    newsp1=sp1                & newsp2=sp2
    newsp1[fcol,*] = newflux1 & newsp2[fcol,*] = newflux2
    newsp1[ecol,*] = newerr1  & newsp2[ecol,*] = newerr2

;   correct the spectra for spikes with calls to spspike.pro

    newsp1=spspike(newsp1,newsp2,_extra=e)
    newsp2=spspike(newsp2,newsp1,_extra=e)

;   combine spectral arrays and sort by wavelength

    testlam  = [lam1,     lam2]
    newsp    = [[newsp1],[newsp2]]
    ordersp  = newsp

    sdx=sort(testlam)

;   load flux and error columns

    ordersp = newsp[*,sdx]

  endif else begin                      ; nothing happens if neither valid
    if (mdx1[0] gt -1) then ordersp=sp1 ; only sp1 valid
    if (mdx2[0] gt -1) then ordersp=sp2 ; only sp2 valid
  endelse

; load outsp array, which will be passed back
; note that ordersp must be defined in first pass

  if (mdx1[0] gt -1 or mdx2[0] gt -1) then begin
    if (m eq minorder) then outsp=ordersp else outsp=[[outsp],[ordersp]]
  endif
 
endfor

RETURN,outsp
END
