FUNCTION sporder,sp,order

;  3 Jan 06 created
;
; sporder returns the specified order in a spectrum
;
; INPUT
; sp    - the spectral data array
; order - the order to be saved; others are discarded
; OUTPUT - returns the one order

ocol=3

idx=where(sp[ocol,*] eq order)
if (max(idx) gt -1) then spout=sp[*,idx] else begin
  print,'Warning, no data in requested order, returning full spectrum.'
  spout=sp
endelse

RETURN,spout
END
