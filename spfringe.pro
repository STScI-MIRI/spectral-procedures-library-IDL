FUNCTION spfringe,lam,amplitude,period,reflam,phase

;  3 Jan 06 created
;
; spfringe generates fringes using an expanding cosine wave
; a baseline of one is assumed, since the fringing is a transmission function
; a cosine is used to assign a phase of zero to a maximum
;
; INPUT
;   lam        - simple vector containing the wavelength grid
;   amplitude  - the amplitude of the fringes about a baseline of 1.0
;   period     - the period of the fringes, in units of the input lam vector
;   reflam     - the wavelength where period and amplitude = first coefficient
;                and the phase=0 (unless set differently)
;   phase      - the phase at reflam (if set) or zero 
;                (runs from 0 to 1, not 2 PI)
; OUTPUT       - a spectral array containing lam and the fringe pattern

lcol=0 & fcol=1 & ncol=4

; update input values for amplitude and period if one-dimensional

if (n_elements(amplitude) eq 1) then amplitude=[amplitude,0]
if (n_elements(period)    eq 1) then period   =[period   ,0]
if (n_elements(reflam)    eq 0) then reflam=0.0
if (n_elements(phase)     eq 0) then phase=0.0

; generate fringes
; first create argument of cosine function
; then modify it to force phase to set value at the reference wavelength reflam

lam=reform(lam)
amp=poly(lam-reflam,amplitude)
per=poly(lam-reflam,period)
argument = (2*!PI/per) * (lam-reflam) + phase

fringe=1.0+amp*cos(argument)

; generate spectral array

spout=fltarr(ncol,n_elements(lam))
spout[lcol,*]=lam
spout[fcol,*]=fringe

RETURN,spout
END
