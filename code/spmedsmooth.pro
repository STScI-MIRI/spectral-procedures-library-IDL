FUNCTION spmedsmooth,sp_in,smoothbox,ORDER=order,NOORDER=noorder,_extra=e

; 19 Jun 12 return input spectrum if smootbox=1
; 12 Jul 10 renamed spmedsmooth
; 21 Dec 06 copied from spsmooth and modified
; 26 Dec 03 created as spsmooth
;
; spmedian breaks a spectrum into orders, median smooths them with the
; IDL median function, reassembles it, and passes it back
;
; the ORDER parameter will limit the smoothing to one order
; if NOORDER keyword set, spectrum will be sorted, smoothed, unsorted
;
; set column identifications, load stops, new sp_array
; assumed that the col 0 = lambda, col 1 = flux, col 2 = error, col 3 = order

lcol=0 & fcol=1 & ecol=2 & ocol=3
ncols=n_elements(sp_in[*,0])
nrows=n_elements(sp_in[0,*])
sp_array  = sp_in

; read keywords
; if order column not present, then set noordflag

if (keyword_set(noorder) eq 1 or nrows lt 3) then noordflag=1 else noordflag=0
if (n_elements(order) eq 0) then order=-1

; find min and max order numbers

minorder=min(sp_array[ocol,*])
maxorder=max(sp_array[ocol,*])
if (order ge 0) then begin                 ; limiting loop to one order
  minorder=order
  maxorder=order
endif

; if noordflag set, then set maxorder=minorder so that main loop runs once
; and sort input array

if (noordflag eq 1) then begin
  maxorder=minorder                        
  idx=sort(reform(sp_in[lcol,*]))          
  sp_array=sp_in[*,idx]
endif
smo_array = sp_array
out_array = sp_array

; loop through orders, median smooth each order and load into smo_array
; error propagation:  smooth square of errors, take sqrt, divide by root N
; error propagation unchanged from spsmooth - still calls smooth, not median

i=0 ; i points to elements in smoothbox, if it's an array
for m=minorder,maxorder do begin

  if (n_elements(smoothbox) gt 1) then smooval=smoothbox[i] $
    else smooval=smoothbox

  if (noordflag eq 0) then ordidx = where(sp_array[ocol,*] eq m) $
    else                   ordidx = indgen(nrows) ; ordidx=everything

  if (max(ordidx) gt -1) then begin
    pass_spec = reform(sp_array[fcol,ordidx])
    err_spec  = sp_array[ecol,ordidx]^2
    if (smooval gt 1) then begin
      new_spec  = median(pass_spec,smooval,_extra=e,/even)
      new_error = sqrt(smooth(err_spec,smooval,_extra=e,/nan)/float(smooval))
    endif else begin
      new_spec = pass_spec
      new_error = err_spec
    endelse
    smo_array[fcol,ordidx] = new_spec
    smo_array[ecol,ordidx] = new_error
  endif

  i=i+1
endfor

; if array sorted, then resort and return, else just return

if (noordflag eq 1) then out_array[*,idx]=smo_array else out_array=smo_array
RETURN,out_array
END
