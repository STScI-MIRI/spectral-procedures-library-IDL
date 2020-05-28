FUNCTION splamshift,spin,dellam,ORDER=order

; 28 Jun 13 created
;
; shifts the wavelengths in spin by dellam (assumed same units)
; if order invoked, only the one order is shifted - others returned unchanged

lcol=0 & fcol=1 & ecol=2 & ocol=3
spout=spin
len=n_elements(spout)

if (keyword_set(order) eq 1) then idx=where(spout[ocol,*] eq order) $
  else idx=indgen(len)

if (idx[0] ge 0) then spout[lcol,idx] += dellam

RETURN,spout
END
