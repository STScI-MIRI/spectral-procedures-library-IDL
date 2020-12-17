FUNCTION splsm,sp,SIGMA=sigma

; 21 Apr 06 created
;
; splsm uses a cubic spline to smooth a spectrum
; it sorts the spectrum and iterates through it
; it replaces each wavelength element with the spline-fit to the rest
; it then sorts the spectrum back to its original order and returns it

lcol=0 & fcol=1 & ecol=2
if (n_elements(sigma) eq 0) then sigma=1.0

len=n_elements(sp[lcol,*])

; sort into array spin

idx=sort(sp[lcol,*])
spin=sp[*,idx]
spout=spin
newlam =reform(spin[lcol,*]) ; will use newlam in spline calls

; iterate through wavelength elements

for i=0,len-1 do begin

; load lam, flux, err arrays (without active element)

  case i of 
    0 : begin
      lam =spin[lcol,1:len-1]
      flux=spin[fcol,1:len-1]
      err =spin[ecol,1:len-1]
    end
    len-1 : begin
      lam =spin[lcol,0:len-2]
      flux=spin[fcol,0:len-2]
      err =spin[ecol,0:len-2]
    end
    else : begin
      lam =[reform(spin[lcol,0:i-1]),reform(spin[lcol,i+1:len-1])]
      flux=[reform(spin[fcol,0:i-1]),reform(spin[fcol,i+1:len-1])]
      err =[reform(spin[ecol,0:i-1]),reform(spin[ecol,i+1:len-1])]
    end
  endcase

; call spline, reload flux and error

  newflux=spline(lam,flux,newlam,sigma)
  newerr =spline(lam,err ,newlam,sigma)
  spout[fcol,i]=newflux[i]
  spout[ecol,i]=newerr[i]

endfor

; unsort spout into spin and return

spin[*,idx]=spout
RETURN,spin
END
