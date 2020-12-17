FUNCTION sptfit,sp,l0,l1,l2,l3,ENGELKE=engelke,DIAG=diag

; 22 May 15 modified to prevent T < 0 in a new iteration
;  5 Sep 07 added diagnostic keyword and improved algorithm
; 31 Aug 07 added Engelke keyword
; 18 Aug 07 created
;
; sptfit fits a Planck or Engelke f'n to a spectrum and returns the temperature
; fit is made over the intervals l0-l1, l2-l3 using spcolor
; the engelke keyword changes the fitting function to an Engelke f'n
;
; INPUT
;   sp      - spectral data array
;   l0,l1   - first fitting interval
;   l2,l3   - second fitting interval
;   engelke - optional keyword to fit Engelke f'n instead of a Planck f'n
;   diag    - turns on diagnostics
;
; OUTPUT - returns temperature
;
; algorithm
;   first find color of spectrum
;   then iterate through temperatures to find best fit to color
;   iteration is by magnitudes, first in intervals of 10,000
;     then in intervals of 1,000 from min +/- 10,000
;     then in 100s, etc.

if (keyword_set(engelke) eq 0) then eflag=0 else eflag=1

; find color of supplied spectrum

color=spcolor(sp,l0,l1,l2,l3)

if (keyword_set(diag) eq 1) then print,color

; set up for iteration

interval=[10000,1000,100,10,1,0.1]

; iterate through intervals

for i=0,n_elements(interval)-1 do begin

; set up stops for temperature range

  if (i eq 0) then begin
    tstart=100
    tstop=1e6
  endif else begin
    tstart = mint - interval[i-1]
    tstop  = mint + interval[i-1]
  endelse

  if (tstart lt 0) then tstart=interval[i]/10.0

; set up t and diff arrays 

  len=(tstop-tstart)/float(interval[i]) +1
  t = tstart + interval[i]*findgen(len)
  diff=dblarr(len)

; iterate through temperatures and fill diff array

  for j=0,len-1 do begin
    if (eflag eq 0) then bb=spplanck(sp,t[j]) else bb=spengelke(sp,t[j])
    testcolor=spcolor(bb,l0,l1,l2,l3)
    diff[j]=abs(testcolor-color)
  endfor

; find minimum diff and corresponding temperature

  mindiff=min(diff,imin)
  mint=t[imin]

; plot and print if in diagnostic mode

  if (keyword_set(diag) eq 1) then begin
    plot,t,diff
    print,i,mint,mindiff
    zq=get_kbrd(1)
  endif

endfor

RETURN,mint
END
