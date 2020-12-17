FUNCTION spnorm,spin,norms

; 29 Dec 04 created
;
; spnorm will multiply each order of spin by the corresponding scalar in norms
;
; for the moment, expecting a one-to-one correspondence with order number
;
; INPUT
;   spin - input spectral data array
;   norms - array of scalars, one per order, plus leading zeroes if needed
;           e.g., if spin is SL, then norms=[0,norm_SL1,norm_SL2,norm_SL-bonus]
; OUTPUT - returns corrected spectrum

ocol=3
spout=spin

minorder=min(spin[ocol,*])
maxorder=max(spin[ocol,*])

for m=minorder,maxorder do begin
  idx=where(spin[ocol,*] eq m)
  if (max(idx) gt -1) then spout=sptimes(spout,norms[m],order=m)
endfor

RETURN,spout
END
