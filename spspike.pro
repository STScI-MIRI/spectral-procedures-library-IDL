FUNCTION spspike,insp1,insp2,SIGMA=sigma,WIDTH=width
;
; 14 Dec 11 created
;
; identify and remove spikes from the flux column of a spectral array
; needs another array on the same or similar grid for comparison
; uses the algorithm developed for sppair.pro (cf for details)
; NOTE - input spectrum must have only one order
;
; INPUT
;   insp1 - input spectrum
;   sigma - keyword to set threshold for identifying spikes and divots
;   width - keyword to set smoothing box size for call to median
; OUTPUT
;   returns the cleaned spectrum

; setup

lcol=0 & fcol=1 & ecol=2 & ocol=3

; check keywords

if (keyword_set(sigma) eq 0) then sigma=5.0
if (keyword_set(width) eq 0) then width=5

; load useful vectors and values

sp1   = insp1
sp2   = spregrid(insp2,insp1)
len   = n_elements(reform(sp1[lcol,*]))
flux1 = reform(sp1[fcol,*])
flux2 = reform(sp2[fcol,*])
con1  = median(flux1,width)
con2  = median(flux2,width)
diff1 = abs(flux1-con1)
diff2 = abs(flux2-con2)
diffx = 0.5*abs(con1-con2)
pass  = fltarr(len)

; generate pass array, min of diff1,diff2 for each pixel

for i=0,len-1 do pass[i] = min([diff1[i],diff2[i]])

; raise low values of pass to median of pass array 

minpass = median(pass)
bad_idx = where(pass lt minpass)
if (bad_idx[0] gt -1) then pass[bad_idx] = minpass

; identify bad pixels if both of the following conditions fulfilled
;   1) pixel is a spike (diff > sigma*pass)
;   2) pixel deviates more than the one in the other nod

bad_idx = where(diff1 gt sigma*pass and diff1 gt diffx)

; replace bad pixels with values from median-smoothed spectrum

if (bad_idx[0] gt -1) then sp1[fcol,bad_idx] = con1[bad_idx]

; return spike-corrected array

RETURN,sp1
END
