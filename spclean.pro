FUNCTION spclean,sp_array,_extra=e

; 26 Dec 03 created
;
; spclean passes a spectrum, order by order to fix_spikes.pro
; in order to clean the spectrum of spikes in the data
; it then reassembles the spectrum in its original order
; and passes it back
;
; assumed that the fourth column has index information

out_array=sp_array

; find min and max order numbers

minorder=min(sp_array[3,*])
maxorder=max(sp_array[3,*])

; loop through orders, clean each order and reload into out_array

for m=minorder,maxorder do begin
  ordidx=where(sp_array[3,*] eq m)
  if (max(ordidx) gt -1) then begin
    pass_spec=sp_array[1,ordidx]
    out_array[1,ordidx]=fix_spikes(pass_spec,_extra=e)
  endif
endfor

return,out_array
END
