PRO spoverplot,inlist,_extra=e,cr=CR

; 28 Feb 08 created
;
; given a file with a list of spectral FITS files
; plot the first and overplot the rest in multiple colors
; this procedure is tailored for a color scheme using colors 2-14

openr,fi,inlist,/get_lun
line=' '

count=0
while (not eof(fi)) do begin

  readf,fi,line
  sp=readfits(line,/silent)
  if (count eq 0) then begin
    spplot,sp,_extra=e 
    count=2
  endif else begin
    spplot,sp,_extra=e,/over,col=count
    if (count le 13) then count=count+1 else count=2
  endelse

  if (keyword_set(cr) ne 0) then zq=get_kbrd(1)  
endwhile

close,/all
END
