FUNCTION spgauss,spec,ga,_extra=e

; 19 Apr 05 added support for nterms=4 and nterms=5
; 15 Apr 05 created
; 
; spgauss creates parameters for 2 gaussians and overplots them
; on a plot of the spectrum "spec"

lcol=0 & fcol=1

nlam=n_elements(spec[lcol,*])
gaussa=spec
for i=0,nlam-1 do $
  gaussa[fcol,i]=ga[0]*exp(-((gaussa[lcol,i]-ga[1])^2.0)/(2.0*ga[2]^2.0))

newfloor=spec 
newfloor[fcol,*]=0.0
if (n_elements(ga) eq 4) then newfloor[fcol,*]=ga[3]
if (n_elements(ga) eq 5) then newfloor[fcol,*]=ga[3]+newfloor[lcol,*]*ga[4]
  
RETURN,spadd(gaussa,newfloor)
END
