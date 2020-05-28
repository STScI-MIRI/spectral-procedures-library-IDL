PRO spgplot,spec,ga,_extra=e
;
; 15 Apr 05
; 
; sp2gplot takes gaussian parameters for 2 gaussians and overplots them
; on a plot of the spectrum "spec"

lcol=0 & fcol=1

nlam=n_elements(spec[lcol,*])
gaussa=spec
for i=0,nlam-1 do gaussa[fcol,i]=gaussg(gaussa[lcol,i],ga[1],ga[2],ga[0])

spplot,spec,_extra=e
spplot,gaussa,/over,th=2,_extra=e

END
