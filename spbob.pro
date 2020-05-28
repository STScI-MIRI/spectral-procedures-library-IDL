FUNCTION spbob,spin

;  2 Mar 12 created
;
; removes 1 pixel from each end of each order of input spectrum
; if the orders are not sorted by wavelength, these pixels won't be real ends
; 
; INPUT
;   spin - input spectrum
; OUTPUT - returns bobbed spectrum

; define columns

ocol=3

omin=min(spin[ocol,*])
omax=max(spin[ocol,*])

count=0
for m=omin,omax do begin

  idx=where(spin[ocol,*] eq m)

  if (idx[0] gt -1) then begin
    sp1=spin[*,idx]
    len=n_elements(idx)
    sp2=sp1[*,1:len-2]
    if (count eq 0) then spout=sp2 else spout=[[spout],[sp2]]
    count=count+1
  endif

endfor

RETURN,spout
END
