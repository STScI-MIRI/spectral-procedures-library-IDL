FUNCTION spsort,spin,NODUPE=nodupe,SINGLEORDER=singleorder

; 20 Dec 12 added singleorder keyword to replace all order values with minimum
; 25 Nov 12 changed loop counters to long to handle big files
; 17 Oct 05 created
;
; sorts a spectral data array by wavelength
; if nodupe keyword set, then average duplicate wavelength elements
;   (this is really the meat of the routine)
;
; INPUT
; spin        - input spectrum
; nodupe      - optional keyword to average all data at the same the
;               wavelength into one wavelength element
; singleroder - if set, will replace all order values with minimum in input

lcol=0 & fcol=1 & ecol=2 & ocol=3

idx=sort(spin[0,*])
spsorted=spin[*,idx]

idx=where(spsorted[lcol,*] eq shift(spsorted[lcol,*],-1) or $
          spsorted[lcol,*] eq shift(spsorted[lcol,*], 1) )

if (keyword_set(nodupe) eq 1 and max(idx) gt -1) then begin

; average flux and propagate errors for duplicate entries

  len=n_elements(idx)
  for i=long(0),len-1 do begin
    mat_idx=where(spsorted[lcol,*] eq spsorted[lcol,idx[i]])
    spsorted[fcol,idx[i]]=mean(reform(spsorted[fcol,mat_idx]))
    spsorted[ecol,idx[i]]= $
      sqrt(mean(reform(spsorted[ecol,mat_idx])^2)/(n_elements(mat_idx)-1))
  endfor

; remove duplicate entries

  len=n_elements(spsorted[lcol,*])
  oldlam=0.0
  spout=spsorted[*,0]
  for i=long(1),len-1 do begin
    if (spsorted[lcol,i] ne oldlam) then begin
      oldlam=spsorted[lcol,i]
      spout=[[spout],[spsorted[*,i]]]
    endif
  endfor
endif else spout=spsorted 

; if singleorder set, then reset order column to minimum of input

if (keyword_set(singleorder) ne 0) then spout[ocol,*]=min(spout[ocol,*])

RETURN,spout
END
