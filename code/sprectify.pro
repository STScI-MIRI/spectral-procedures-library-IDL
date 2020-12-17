FUNCTION sprectify,spin,lambda,NLAM=nlam,DEGREE=degree,BOX=box,PLOT=plot

; 16 Apr 13 created
;
; sprectify straightens a spectrum with a quadratic
; NOTE - not treating orders separately!
;
; INPUT
;   spin - input spectrum
;   lambda - optional set of wavelength grid points
;   nlam   - keyword to set number of grid points if lambda not supplied
;            DEFAULT=10
;   degree - degree of polynomial to fit, DEFAULT=4
;   box    - size of smoothing box, DEFAULT=3

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check keywords
; construct lambda vector if necessary

if (keyword_set(degree) eq 0) then degree=4
if (keyword_set(box) eq 0) then box=3

if (keyword_set(lambda) eq 0) then begin
  if (keyword_set(nlam) eq 0) then nlam=10
  minlam=min(spin[lcol,*])
  maxlam=max(spin[lcol,*])
  lambda = minlam + findgen(nlam)*(maxlam-minlam)/(nlam-1)
endif

; sort data by wavelength (ignoring orders)

idx=sort(spin[lcol,*])
sp=spin[*,idx]

; smooth the spectrum

if (box gt 1) then smoothed=spsmooth(sp,box)

; fit a polynomial to the data at the lambda positions

flux=lambda
for i=0,nlam-1 do flux[i] = spspot(smoothed,lambda[i])
cc=poly_fit(lambda,flux,degree)

; rectify the spectrum with the polynomial - in array spc

spc=sp
correction=poly(sp[lcol,*],cc)
spc[fcol,*]=sp[fcol,*]/correction
spc[ecol,*]=sp[ecol,*]/correction

; plot if requested

if (keyword_set(plot) eq 1) then begin
  spplot,sp
  oplot,lambda,flux,psym=4,symsize=2
  oplot,sp[lcol,*],correction
endif

; resort to return to original order

spout=spin
spout[*,idx]=spc

RETURN,spout
END
