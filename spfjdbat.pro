PRO spfjdbat,inlist,outfile,_extra=e

; 30 Jul 11 created
;
; reads inlist, read and prints the MJD_OBS keyword from the FITS header
; expects one file per line in inlist

openr,fi,inlist,/get_lun
line=' '
fmt='(a15,1x,f12.5)'

if (n_elements(outfile) gt 0) then outflag=1 else outflag=0
if (outflag eq 1) then openw,fo,outfile,/get_lun

while (not eof(fi)) do begin

  readf,fi,line

  split1=strsplit(line,' ',/extract)
  split2=strsplit(split1[0],'.',/extract)
  target=split2[0]+'.'+split2[1]
  while (strlen(target) lt 15) do target=target+' '

  data=readfits(split1[0],hdr,/silent)
  mjd=sxpar(hdr,'MJD_OBS')

  if (outflag eq 1) then printf,fo,target,mjd,format=fmt else $
    print,target,mjd,format=fmt

endwhile

free_lun,fi
if (outflag eq 1) then free_lun,fo
END
