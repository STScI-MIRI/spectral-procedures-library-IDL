FUNCTION spgridres,spold,spnew,fwhm,RANGE=range,VERBOSE=verbose

;  8 Jan 14 debugging
; 30 Dec 13 rewritten as spgridres with FWHM a polynomial f'n of wavelength
;           FWHM is in units of um, vs. the old sigma in pixel units
;           range is still input as a multiple of sigma (!)
;  1 Nov 12 added verbose keyword, modified to handle 1-D input for spnew
; 30 Jun 05 modified from gridlo.pro, version 4.0
; need to include error analysis, but that will come later
;
; gridlo.pro history
;  3 Oct 03 modify range to be in units of sigma and sigma to be in 
;           units of pixel (in the old or new wavelength grid, as set by OLD)
;           also take delta lambda to be larger of result to either side
; 10 Sep 03 add keyword old, if set, use old wavelength grid to find sigma
;           also corrected a bug which forced sigma to a constant value
; 31 Jan 02 change so sigma set based on new grid, not old 
; 24 Jan 02 tidy up code and bugs
; 11 Dec 01 simplify Russ's modifications
; 31 Oct 01 Russ to speed it up with where commands
;  1 Oct 01 make counters long
; 27 Mar 97
; 29 May 93 created
;
;  this procedure converts spectra from one wavelength grid to another
;  convolves a gaussian centered over each new
;    pixel to weight data from old spectrum
;  old warning message:
;    IF YOU'RE GOING FROM LOW to HIGH RES, use REGRID (or gridhi) instead!
;    this was because the old lo-res grid would leave data gaps in the new
;    with the Dec 01 - Jan 02 improvements, there should be no gaps
;  outstanding issue, what to use for sigma if not declared
;    changed 31 Jan 02 to set sigma based on new wavelength grid
;
;  input
;    spold   - spectrum to be regridded
;              a = input spectrum, need at least 2 columns, wavelength and flux
;    spnew   - spectrum with new wavelength grid - can now be 1-D
;              l = new wavelength grid, loaded from spnew
;    range   - number of sigma spacing to integrate over (default=3)
;              units are multiples of sigma
;    fwhm    - FWHM of gaussian to integrate over (in um)
;              this is input as a polynomial f'n of wavelength
;    verbose - if set, then each wavelength printed to screen as it's finished
;  output   flux on new wavelength grid

; initialize 

lcol=0 & fcol=1 & ecol=2 & ocol=3
flipflag=0

; read in keywords and set defaults

if (keyword_set(fwhm) eq 0) then fwhmflag=0  else fwhmflag=1
if (keyword_set(range) eq 0) then range=3.0
if (keyword_set(verbose) eq 0) then verbflag=0 else verbflag=1

; sort spold into wavelength order

data = spold
idx  = sort(data[0,*])
data = spold[*,idx]

; load l vector = new wavelength grid (unsorted)

sz=size(spnew)
case sz[0] of
  1 : begin & ncols=1                      & newlam=spnew                 & end
  2 : begin & ncols=n_elements(spnew[*,0]) & newlam=reform(spnew[lcol,*]) & end
  else : print,'Crash coming - need 1-D or 2-D input for new grid'
endcase

; setup continued - now just using old gridlo code

oldlam  = reform(data[0,*])
oldflux = reform(data[1,*])
len  = long(n_elements(oldlam))
oldstop = len-1
newstop = n_elements(newlam)-1
newflux = fltarr(newstop+1)

; reverse order of al and af if al is not a positive grid
; ditto for l, the new wavelength grid

if (oldlam[0] gt oldlam[oldstop]) then begin
  oldlam  = rotate(oldlam,2)
  oldflux = rotate(oldflux,2)
endif
if (newlam[0] gt newlam[newstop]) then begin
  newlam=rotate(newlam,2)
  flipflag=1 
endif

; set grid of deltas for new wavelength grid (newdel) and old grid (olddel) 
; check delta to either side, take larger

newdel=newlam
for j=1,long(newstop)-1 do $
  newdel[j]     = max([newlam[j+1]-newlam[j],newlam[j]-newlam[j-1]])
newdel[0]       = newlam[1]-newlam[0]
newdel[newstop] = newlam[newstop]-newlam[newstop-1]

olddel=oldlam
for j=long(1),oldstop-1 do $
  olddel[j]     = max([oldlam[j+1]-oldlam[j],oldlam[j]-oldlam[j-1]])
olddel[0]       = oldlam[1]-oldlam[0]
olddel[oldstop] = oldlam[oldstop]-oldlam[oldstop-1]

; if FWHM not defined, then define it so that sigma = spacing in new grid
; generate a sigma vector from FWHM on the new wavelength grid

if (fwhmflag eq 0) then sigma=newdel/2.354 else sigma=poly(newlam,fwhm)/2.354

; iterate through wavelength elements of new grid
; old and new wavelength grids are now positive

for i=long(0),newstop do begin

  sum=0.0
  wtsum=0.0

; the algorithm is
;   (1) find the position in the old grid corresponding to the new position
;   (2) set sigma and integration range about this point
;   (3) perform the weighted gaussian summation about this point
;
; (1) find alspot (set to 0 if new grid point too low)

  lam     = newlam[i]
  oldspot = max(where(lam gt oldlam)) > 0

; (2) find sig and ran - relevant values of sigma and range*sigma

  sig = sigma[i]
  ran = range*sig

; (3) find range of old grid counters within range of alspot
; oldrange gives the indices over which to integrate

  oldrange = where(oldlam ge lam-ran and oldlam le lam+ran,oldcount)
  minrange = min(oldrange)
  maxrange = max(oldrange)

; (4) if oldcount gives 0 pixels in range, pick oldspot 
; and neighbors if oldspot not at the end of the range

  if (oldcount eq 0) then begin
    oldrange = oldspot
    oldcount = 1
    if (oldspot lt oldstop) then begin
      oldrange = [oldrange, oldspot+1]
      oldcount ++
    endif 
    if (oldspot gt 0) then begin
      oldrange = [oldspot-1, oldrange]
      oldcount ++
    endif
  endif

; (5) add one data point to each end of range, unless at 0 or last point
; this would guarantee at least 3 data points at ends, and 5 in middle

  if (oldcount eq 0) then print,'error in spgridres at ',i ; shouldn't happen
  if (minrange gt 0 and maxrange lt oldstop) then begin
    oldrange = [minrange-1,oldrange,maxrange+1]
    oldcount += 2
  endif else begin
    if (minrange eq 0) then begin
      oldrange = [oldrange,maxrange+1]
      oldcount ++
    endif
    if (maxrange eq oldstop) then begin
      oldrange = [minrange-1,oldrange]
      oldcount ++
    endif
  endelse

; (6) integrate:  gaussian weighted sum over range

  if (oldcount gt 0) then begin
    wt    = exp(-((oldlam[oldrange]-lam)*(oldlam[oldrange]-lam))/(2.0*sig*sig))
    sum   = total(wt*oldflux[oldrange])
    wtsum = total(wt)

    if (wtsum gt 0.0) then newflux[i] = sum/wtsum else newflux[i] = 0.0

;   load flux array

  endif else begin
    newflux[i] = 0.0
    print,'no data in range at pixel ',i
  endelse
  if (newflux[i] eq 0.0 and i lt 1000 and verbflag eq 1) then $
    print,'zero at pixel ',i

; (7) say something if /verbose set

  if (verbflag eq 1) then print,i,lam,format='(i8,f15.5)'
endfor

; load newflux as second column of returned spectral array

case ncols of
  1 : begin
      spreturn=fltarr(4,n_elements(newlam))
      spreturn[lcol,*]=newlam               ; load wavelength
      end
  2 : spreturn=spnew                        ; automatically loads wavelength
  else : begin
      spreturn=spnew
      spreturn[ecol,*]=0.0                  ; zero errors out
      end
endcase

; check to see if we flipped wavelength grid - if so, flip newflux back

if (flipflag eq 1) then newflux=rotate(newflux,2)

spreturn[fcol,*]=newflux                    ; load new spectrum 

RETURN,spreturn
END
