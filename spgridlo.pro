FUNCTION spgridlo,spold,spnew,RANGE=range,SIGMA=sigma,OLD=old,VERBOSE=verbose

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
;    sigma   - sigma of gaussian to integrate over (in lambda units)
;              units are multiples of wavelength spacing between pixels
;    old     = 1 to use old grid to determine sigma, 0 for new grid (default=0)
;    verbose - if set, then each wavelength printed to screen as it's finished
;  output:  flux on new wavelength grid

; initialize 

lcol=0 & fcol=1 & ecol=2 & ocol=3

; read in keywords and set defaults

if (n_elements(sigma) eq 0) then begin
  sigmaflag=0 
endif else begin
  if (float(sigma) eq 0) then sigmaflag=0 else sigmaflag=1
endelse

if (keyword_set(old) eq 0) then old=0 else old=1
if (keyword_set(range) eq 0) then range=3.0
if (float(range) eq 0.0) then range=3.0
if (keyword_set(verbose) eq 0) then verbflag=0 else verbflag=1

; sort spold into wavelength order

a=spold
idx=sort(a[0,*])
a=spold[*,idx]

; load l vector = new wavelength grid (unsorted)

sz=size(spnew)
case sz[0] of
  1 : begin & ncols=1                      & l=spnew                 & end
  2 : begin & ncols=n_elements(spnew[*,0]) & l=reform(spnew[lcol,*]) & end
  else : print,'Crash coming - need 1-D or 2-D input for new grid'
endcase

; setup continued - now just using old gridlo code

al=reform(a[0,*])
af=reform(a[1,*])
nal=long(n_elements(al))
alstop=nal-1
lstop=n_elements(l)-1
b=fltarr(lstop+1)

; reverse order of al and af if al is not a positive grid
; ditto for l, the new wavelength grid

if (al[0] gt al[alstop]) then begin
  al=rotate(al,2)
  af=rotate(af,2)
endif
if (l[0] gt l[lstop]) then l=rotate(l,2)

; set grid of deltas for new wavelength grid (del) and old grid (dal) 
; check delta to either side, take larger

del=l
for j=1,lstop-1 do del[j]=max([l[j+1]-l[j],l[j]-l[j-1]])
del[0] = l[1]-l[0]
del[lstop] = l[lstop]-l[lstop-1]

dal=al
for j=long(1),alstop-1 do dal[j]=max([al[j+1]-al[j],al[j]-al[j-1]])
dal[0] = al[1]-al[0]
dal[alstop] = al[alstop]-al[alstop-1]

; iterate through wavelength elements of new grid
; old and new wavelength grids are now positive

for i=long(0),lstop do begin

  sum=0.0
  wtsum=0.0

; the algorithm is
;   (1) find the position in the old grid corresponding to the new position
;   (2) set sigma and integration range about this point
;   (3) perform the weighted gaussian summation about this point
;
; (1) find alspot (set to 0 if new grid point too low)

  lam = l[i]
  alspot=max(where(lam gt al)) > 0

; (2a) set sigma and range if not already set (i.e. = 0.0)
; if sigma = 0.0, set so FWHM = wavelength spacing in the old grid

  if (sigmaflag eq 0) then begin
    if (old eq 0) then sig = 0.424 * del[i] $
    else               sig = 0.424 * dal[alspot]
  endif else begin
    if (old eq 0) then sig = sigma * del[i] $
    else               sig = sigma * dal[alspot]
  endelse

  ran = range*sig

; (2b) find range of old grid counters within range of alspot
; alrange gives the indices over which to integrate

  alrange = where(al ge lam-ran and al le lam+ran,alcount)
  minrange=min(alrange)
  maxrange=max(alrange)

; if alcount gives 0 pixels in range, pick alspot 
; and neighbors if alspot not at the end of the range

  if (alcount eq 0) then begin
    alrange = alspot
    alcount = 1
    if (alspot lt alstop) then begin
      alrange = [alrange, alspot+1]
      alcount = alcount+1
    endif 
    if (alspot gt 0) then begin
      alrange = [alspot-1, alrange]
      alcount = alcount+1
    endif
  endif

; add one data point to each end of range, unless at 0 or last point
; this would guarantee at least 3 data points at ends, and 5 in middle

  if (alcount eq 0) then print,'error in gridlo at ',i ; this shouldn't happen
  if (minrange gt 0 and maxrange lt alstop) then begin
    alrange = [minrange-1,alrange,maxrange+1]
    alcount = alcount+2
  endif else begin
    if (minrange eq 0) then begin
      alrange = [alrange,maxrange+1]
      alcount = alcount+1
    endif
    if (maxrange eq alstop) then begin
      alrange = [minrange-1,alrange]
      alcount = alcount+1
    endif
  endelse

; (3) integrate:  gaussian weighted sum over range

  if (alcount gt 0) then begin
    wt = exp(-((al[alrange]-lam)*(al[alrange]-lam))/(2.0*sig*sig))
    sum=total(wt*af[alrange])
    wtsum=total(wt)

    if (wtsum gt 0.0) then b[i] = sum/wtsum else b[i] = 0.0

;   load flux array

  endif else begin
    b[i] = 0.0
    print,'no data in range at pixel ',i
  endelse
  if (b[i] eq 0.0 and i lt 1000 and verbflag eq 1) then print,'zero at pixel ',i

; (4) say something if /verbose set

  if (verbflag eq 1) then print,i,lam,format='(i8,f15.5)'
endfor

; load b as second column of returned spectral array

case ncols of
  1 : begin
      spreturn=fltarr(4,n_elements(l))
      spreturn[lcol,*]=l   ; load wavelength
      end
  2 : spreturn=spnew       ; automatically loads wavelength
  else : begin
      spreturn=spnew
      spreturn[ecol,*]=0.0 ; zero errors out
      end
endcase

spreturn[fcol,*]=b         ; load new spectrum 

RETURN,spreturn
END
