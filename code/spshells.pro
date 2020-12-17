FUNCTION spshells,q,t,f

; 28 Apr 06 modified so that "t" is an array of temperatures
; 26 Apr 06 created
;
; spshells generates a series of dust shells and returns the spectrum
; the optically thin case is assumed
;
; INPUT
;   q    - array with 2 or more columns, this defines the dust
;          col 0 = wavelength (um), col 1 = Q (optical efficiency)
;          this defines the dust
;   t    - an array of temperatures, one per shell
;   f    - an array of fill factors, one for each shell
;          NOTE - negative values are treated as positive!
;
; OUTPUT - spectral array which matches q in dimensions
;          flux column = S(lambda) = sum (f(T) * B(T) * Q(lambda)
;

lcol=0 & fcol=1 & ecol=2

; set up output spectral array
; four columns - wavelength (um), flux (Jy), error (Jy), order
;   wavelength taken from q[0,*] and assumed to be in um
;   flux computed in F_nu space
;   error and order set to zero

nlam  = n_elements(q[lcol,*])
ncol  = n_elements(q[*,0])
nshell= n_elements(t)

spout = q
spout[fcol,*]=0.0
if (ncol gt ecol) then spout[ecol,*]=0.0

for i=0,nshell-1 do begin
  sptemp=spplanck(q,t[i])
  sptemp=sptimes(sptemp,q)
  spout = spadd(spout,sptimes(sptemp,abs(f[i])))
endfor

RETURN,spout
END
