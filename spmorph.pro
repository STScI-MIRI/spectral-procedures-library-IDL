FUNCTION spmorph,in,fit,grid

; 26 Jul 18 created
;
; morphs one spectrum to another using a spline
; find the ratio of fit/in at each gridpoint using spspot (for now)
; fits a spline to the ratios
; output spectrum is on the grid of the input spectrum

; test procedure, may form the basis of more sophisticated code later

; initializations

lcol=0 & fcol=1 & ecol=2

; set up ratio vector

len=n_elements(grid)
ratio = fltarr(len)

; for each element in the ratio vector, find the ratio with spspot calls

for i=0,len-1 do ratio[i] = spspot(fit,grid[i])/spspot(in,grid[i])

; pad the grid and ratio vectors to cover the wavelength range of in spectrum
; all ratios to blue = blue-most ratio, all ratios to red = red-most ratio

; blue first

nblue = fix(min(fit[lcol,*])-min(in[lcol,*]))
if (nblue gt 0) then begin
  bluegrid = fix(min(in[lcol,*])) + indgen(nblue)
  blueratio = replicate(ratio[0],n_elements(bluegrid))
  grid = [bluegrid,grid]
  ratio = [blueratio,ratio]
endif

; now red

nred = fix(max(in[lcol,*])-max(fit[lcol,*]))
if (nred gt 0) then begin
  redgrid = fix(max(in[lcol,*])+1) + indgen(nred)
  redratio = replicate(ratio[len-1],n_elements(redgrid))
  grid = [grid,redgrid]
  ratio = [ratio,redratio]
endif
 
; fit a spline to generate corrections
; may need to check for monotonic grid here

corr = cspline(grid,ratio,reform(in[0,*]))

; set up output spectrum

out = in

; multiply flux density and uncertainty columns by corr vector

out[fcol,*] = in[fcol,*]*corr
out[ecol,*] = in[ecol,*]*corr

; return output spectrum

RETURN,out
END
