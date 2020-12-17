FUNCTION spzap,sparray,l0,l1,order=ORDER

;  2 Nov 12 modified for only 2 columns in sparray 
; 31 Oct 12 modified so that l0 and l1 aren't necessary if order specified
;  9 Jun 06 repaired bug when all of first order zapped
; 20 Jul 04 copied from spcut.pro
;
; spzap removes data in a spectral data array between wavelengths l0 and l1
; INPUT
;   sparray - input spectral data array
;   l0      - minimum wavelength of range to zap
;   l1      - maximum wavelength of range to zap
;   order   - order to zap
;             if set, expects order info in col 3
;             if out of range, then nothing trimmed

; initialize

ncol=n_elements(sparray[*,0])

lcol=0 & ocol=3
if (ncol lt ocol) then ocol=1 ; reset ocol if only 2 cols present

if (n_elements(order) eq 0) then ordflag=0 else ordflag=1
full_len=n_elements(reform(sparray[0,*]))
zapflag=1
emptyflag=0 ; set to one if an order is completely zapped

; check wavelengths, set to zero and 9999.0 if not set but order is

if (ordflag eq 1) then begin
  if (keyword_set(l1) eq 0) then l1=9999.0
  if (keyword_set(l0) eq 0) then l0=0.0
endif

; isolate relevant order as ordarray

if (ordflag eq 1) then begin
  ord_idx=where(sparray[ocol,*] eq order)
  if (ord_idx[0] gt -1) then begin
    ordarray=sparray[*,ord_idx] 
  endif else begin
     print,'Warning in spzap.  No data in order ',$
       order,'.  Returning input array.'
     zapflag=0
  endelse
endif else ordarray=sparray ; no order set, load full input array

; zap data in ordarray, put into savearray

if (zapflag eq 1) then begin
  lam=reform(ordarray[lcol,*])
  out_idx=where( (lam lt l0) or (lam gt l1) )
  if (out_idx[0] gt -1) then savearray=ordarray[*,out_idx] $
  else begin
    emptyflag=1
    savearray=ordarray
  endelse 
endif

; build output array

if (zapflag eq 0) then outarray=sparray ; load full input array
if (zapflag eq 1) then begin

; load orders before active order, active order (as cut), then orders after

  if (ordflag eq 0) then outarray=savearray else begin
    minact=min(ord_idx) & maxact=max(ord_idx)
    if (minact gt 0) then begin
      outarray=sparray[*,0:minact-1]
      if (emptyflag ne 1) then outarray=[[outarray],[savearray]]
      if (maxact lt full_len-1) then $
        outarray=[[outarray],[sparray[*,maxact+1:full_len-1]]]
    endif else begin
      if (emptyflag ne 1) then begin
        outarray=savearray
        outarray=[[outarray],[sparray[*,maxact+1:full_len-1]]]
      endif else outarray=sparray[*,maxact+1:full_len-1]
    endelse
;      if (emptyflag eq 0) then $
;        else outarray=sparray[*,maxact+1:full_len-1]
  endelse

endif

return,outarray

END
