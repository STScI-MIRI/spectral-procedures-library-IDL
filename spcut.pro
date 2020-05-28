FUNCTION spcut,sparray,l0,l1,ORDER=order,ALL=all

;  1 Jan 14 if l0 and l1 not passed, allflag=1
; 15 Jun 06 rewrote last block to properly load various pieces of output array
; 28 Jan 06 corrected error resulting in change to l0 and l1 
; 21 Aug 05 added ALL keyword to cut entire order, must set ORDER keyword, too
; 27 Mar 04 added option to cut only one order at a time
;  1 Feb 04 return error message and orig. data if no data in range
; 26 Jan 04 created
;
; spcut trims a spectral data array to the wavelength range l0-l1
; it is best to apply this to only one spectral segment at a time
; INPUT
;   sparray - input spectral data array
;   l0      - minimum wavelength to save
;   l1      - maximum wavelength to save
;   order   - optional parameter, order to trim
;             if set, expects order info in col 3
;             if out of range, then nothing trimmed
;   all     - optional keyword, works only if order set
;             will trim entire order

; read keywords

lcol=0 & ocol=3
full_len=n_elements(reform(sparray[lcol,*]))
if (n_elements(order) eq 0) then ordflag=0 else ordflag=1
allflag=0
if (ordflag eq 1) then begin
  if (keyword_set(l0) eq 0 and keyword_set(l1) eq 0) then allflag=1
  if (keyword_set(all) eq 1) then allflag=1
endif
cutflag=1     ; set to zero if no cutting to be done
emptyflag=0   ; set to one if active order to be completely removed

; load lambda stops from input l0 and l1 (so that l0 and l1 don't change)

if (allflag eq 0) then begin
  la=l0 & lb=l1
endif else begin
  la=0 & lb=0
endelse

; isolate relevant order as ordarray

if (ordflag eq 1) then begin
  ord_idx=where(sparray[ocol,*] eq order)
  if (max(ord_idx) gt -1) then begin
    ordarray=sparray[*,ord_idx] 
    if (allflag eq 1) then begin ; reset la and lb to max of order +1 and max+2
      la=max(ordarray[lcol,*])+1
      lb=la+1
    endif
  endif else begin
     print,'Warning in spcut.  No data in order ',$
       order,'.  Returning input array.'
     cutflag=0
  endelse
endif else ordarray=sparray ; no order set, load full input array

; perform cut on ordarray, put into savearray

if (cutflag eq 1) then begin
  lam=reform(ordarray[lcol,*])
  out_idx=where( (lam ge la) and (lam le lb) )
  if (max(out_idx gt -1)) then savearray=ordarray[*,out_idx] $
  else savearray=-1       ; -1 value means won't be used below
endif

; build output array

if (cutflag eq 0) then begin
  outarray=sparray                            ; load full input array
endif else begin

  if (ordflag eq 0) then begin                ; load for no set order

    if (max(savearray) gt -1) then outarray=savearray else begin
      print,'Warning in spcut.  No data in range.  Returning input array.'
      outarray=sparray
    endelse

  endif else begin                            ; load when order set

; load orders before active order (firstarray), active order (savearray)
; and orders after (lastarray)

    minact=min(ord_idx) & maxact=max(ord_idx)
    if (minact gt 0) then firstarray=sparray[*,0:minact-1] $
      else firstarray=-1
    if (maxact lt full_len-1) then lastarray=sparray[*,maxact+1:full_len-1] $
      else lastarray=-1

    if (max(firstarray) gt -1) then begin            ; start with firstarray
      outarray=firstarray 
      if (max(savearray) gt -1) then outarray=[[outarray],[savearray]]
      if (max(lastarray) gt -1) then outarray=[[outarray],[lastarray]]
    endif else if (max(savearray) gt -1) then begin ; start with savearray
      outarray=savearray
      if (max(lastarray) gt -1) then outarray=[[outarray],[lastarray]]
    endif else if (max(lastarray) gt -1) then begin ; lastarray only
      outarray=lastarray
    endif else begin                                ; have to load something
      outarray=sparray
      print,'Warning in spcut.  No valid data.  Returning input array.'
    endelse

  endelse

endelse

return,outarray
END
