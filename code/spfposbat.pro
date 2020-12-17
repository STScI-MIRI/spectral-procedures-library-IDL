PRO spfposbat,inlist,outfile,FULL=full,ALL=all,SEX=sex,_extra=e

;  6 Mar 14 added /sex keyword to print RA and DEC in sexegismal units
;           if /sex set, then /all ignored
;           added an extra digit to decimal output
; 15 Feb 13 added /all keyword to return EXTPOS along with RA_EXT and DEC_EXT
; 20 Feb 12 added /full keyword to print full target name (but not path)
; 16 Aug 11 modified to remove path from file name
; 22 Jul 11 created
;
; reads inlist, makes repeated calls to spfpos.pro
; expects one file per line in inlist

openr,fi,inlist,/get_lun
line=' '
if (keyword_set(full) eq 0) then len=15 else len=30
if (keyword_set(sex) eq 0) then degflag=1 else degflag=0
lenstr=string(len,format='(i2)')
fmtd='(a'+lenstr+',1x,2(f11.6),f8.3)'
fmts='(a'+lenstr+',1x,2(a12))'

if (n_elements(outfile) gt 0) then outflag=1 else outflag=0
if (outflag eq 1) then openw,fo,outfile,/get_lun

while (not eof(fi)) do begin

  readf,fi,line

  split1=strsplit(line,' ',/extract)
  split2=strsplit(split1[0],'/',/extract)
  split3=strsplit(split2[n_elements(split2)-1],'.',/extract)
  if (keyword_set(full) eq 0) then target=split3[0]+'.'+split3[1] else $
    target=split2[n_elements(split2)-1]
  while (strlen(target) lt len) do target=target+' '

  data=spfpos(split1[0],_extra=e)
  if (keyword_set(all) ne 0 and degflag eq 1) then $
    data=[data,spfpos(split1[0],/pos)]
  if (outflag eq 1) then begin
    if (degflag eq 1) then printf,fo,target,data,format=fmtd $
      else printf,fo,target,deg2sex(data[0],data[1],/vector),format=fmts
  endif else begin
    if (degflag eq 1) then print,target,data,format=fmtd $
      else print,target,deg2sex(data[0],data[1],/vector),format=fmts
  endelse
endwhile

free_lun,fi
if (outflag eq 1) then free_lun,fo
END
