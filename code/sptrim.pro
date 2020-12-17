FUNCTION sptrim,sp,l0,l1,ORDARRAY=ordarray

;  2 Nov 12 modified to handle 2-column input array
;  9 Dec 11 minor clean-up
;  8 Jun 06 modify to allow a single l0 and l1 to be passed
; 28 Jan 06 modified to cut segments with no data in l0-l1 range
; 17 Sep 05 added ordarray parameter
;  8 Dec 04 modified to be generic
; 28 Dec 03 created
;
; trims bad data from ends of orders
;
; INPUT
;   sp       - spectral data array
;   l0       - vector of starting wavelength for each order
;   l1       - vector of ending wavelength for each order
;   ordarray - vector of order numbers in order as given in l0, l1
;              this parameter helps the program in the event of missing orders
; OUTPUT - returns trimmed spectra data array
;
; assumed that columns are wavelength, flux, error, order

ncol=n_elements(sp[*,0])

lcol=0 & ocol=3
if (ncol lt ocol) then ocol=1

; setup

lam  =reform(sp[lcol,*]) 
order=reform(sp[ocol,*])
m0=min(order) & m1=max(order)

; if ordarray not given, then build it

i=0
if (keyword_set(ordarray) eq 0) then begin
  for m=m0,m1 do begin
    o_idx=where (order eq m)
    if (max(o_idx) gt -1) then begin
      if (i eq 0) then ordarray=m else ordarray=[ordarray,m]
    endif
    i=i+1
  endfor
endif

; expand l0 and l1 to vectors if they are scalars

if (n_elements(l0) eq 1) then begin
 test=l0 & l0=fltarr(n_elements(ordarray)) & l0[*]=test
 test=l1 & l1=fltarr(n_elements(ordarray)) & l1[*]=test
endif

; iterate through orders and trim them
; load into out_array as we go

setflag=0 ; changes to 1 when first order loaded into out_array
for i=0,n_elements(ordarray)-1 do begin
  o_idx=where(order eq ordarray[i])
  if (o_idx[0] gt -1) then begin
    seg_array=sp[*,o_idx]
    olam=lam[o_idx]
    save_idx=where( (olam ge l0[i]) and (olam le l1[i]) )
    if (max(save_idx) gt -1) then begin
      cut_array=seg_array[*,save_idx ]
      if (setflag eq 0) then begin
        out_array=cut_array
        setflag=1
      endif else out_array=[[out_array],[cut_array]]
    endif
  endif else print,'Warning in sptrim, no data in order ',ordarray[i]
endfor

; if setflag still 0, return original array and print error message 

if (setflag eq 0) then begin
  out_array=sp
  print,'Warning in sptrim, returning original spectral array untrimmed.'
endif

return,out_array

END
