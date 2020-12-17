PRO spposstack,infile
;
; 21 Feb 13 created
;
; reads output from spfposbat
; averages RA and DEC for each source
; expects the followign columns:  target, RA, dec

; open input file

openr,fi,infile,/get_lun
line=' '

oldtarget='oldtarget'
fmt='(a20,2(f11.5,f9.5))'

count=0
while (not eof(fi)) do begin

; read and parse one line of input

  readf,fi,line
  if (strlen(line) gt 0) then begin
    split1 = strsplit(line,' ',/extract)
    split2 = strsplit(split1[0],'.',/extract)
    target = split2[0]+'.'+split2[1]
    radum  = float(split1[1])
    decdum = float(split1[2])

;   if new data then start collecting, else report on old data

    if (target ne oldtarget) then begin
      if (oldtarget ne 'oldtarget') then begin          ; report on old data
        while (strlen(oldtarget) lt 20) do oldtarget += ' '
        if (count gt 1) then $
          print,oldtarget,mean(ra),stddev(ra),mean(dec),stddev(dec),format=fmt $
        else print,oldtarget,mean(ra),0.0,mean(dec),0.0,format=fmt 
      endif
      ra=radum & dec=decdum
      oldtarget=target
      count=1
    endif else begin
      ra=[ra,radum] & dec=[dec,decdum]
      count += 1
    endelse
  endif
endwhile

; have to report on last set of data

while (strlen(target) lt 20) do target += ' '
print,target,mean(ra),stddev(ra),mean(dec),stddev(dec),format=fmt 

free_lun,fi
END
