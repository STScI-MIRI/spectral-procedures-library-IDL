PRO dobat,input,counter,count,goflag=goflag,zoom=zoom,$
  xrange=xrange,range=range,filter=filter

; 21 Apr 15 added filter keyword (F and f)
; 17 Nov 11 added '?' for help
; 22 Jan 09 created
;
; dobat reads keyboard input and translates it for a batch display program
; the input typically advances a counter for a list file
;   or updates plotting ranges or similar functions
; INPUT
;   input - keystroke to be read
;   count - number of files in list file
; OUTPUT
;   returns counter value
; MODIFIABLE KEYWORDS
;   counter - position in file list
;   goflag  - flag passed back to main routine, 1 = continue, 2 = conclude
;   range   - 2-element display range (for images, or spectra - as yrange)
;   zoom    - 1-element zoom factor (for images, not spectra)
;   xrr     - 2-element display x range (for spectra, not images)

line=' '

case byte(input) of
    27 : action='quit'   ; ESC
    69 : action='quit'   ; E
    81 : action='quit'   ; Q
   101 : action='quit'   ; e
   113 : action='quit'   ; q
    75 : action='up'     ; K
    85 : action='up'     ; U
   107 : action='up'     ; k
   117 : action='up'     ; u
    74 : action='jump'   ; J
   106 : action='jump'   ; j
    66 : action='up'     ; B
    98 : action='up'     ; b
    84 : action='top'    ; T
   116 : action='top'    ; t
    76 : action='last'   ; L
   108 : action='last'   ; l
    82 : action='range'  ; R
   114 : action='range'  ; r
    90 : action='zoom'   ; Z
   122 : action='zoom'   ; z
    88 : action='xrange' ; X
   120 : action='xrange' ; x
    70 : action='filter' ; F
   102 : action='filter' ; f
    72 : action='help'   ; H
   104 : action='help'   ; h
    63 : action='help'   ; ?
   126 : action='bad'    ; ~
  else : action='down'
endcase

; process action

case action of

  'quit'   : goflag=0

  'up'     : begin
    if (counter gt 0) then counter=counter-1 else begin
      print,'At top of list.  Redisplaying first file.'
    endelse
  end

  'down'   : begin
    if (counter lt count) then counter=counter+1 else begin
      print,"At bottom of list.  Redisplaying last file.  Type 'q' to quit."
    endelse
  end

  'jump'   : begin
    print," "
    print,"Enter a line number, 't' for top, 'l' for last, or +/-decrement"
    read,line
    char=byte(strmid(line,0,1))
    case char of
        76 : counter=count ; L
        84 : counter=0     ; T
       108 : counter=count ; l
       116 : counter=0     ; t
        43 : counter=counter+fix(strmid(line,1)) ; +
        45 : counter=counter-fix(strmid(line,1)) ; -
      else : if (char ge 48 and char le 57) then counter=fix(strmid(line,0))
    endcase
    if (counter lt 0) then counter=0
    if (counter gt count) then counter=count
  end

  'last'   : counter=count

  'top'    : counter=0

  'range'  : begin
   if (n_elements(range) eq 2) then begin
      print,'Display range for signal is ',range[0],' to ',range[1]
   endif else print,'Display range for signal is ', range[0]
    print,'Enter a new (space-delimited) range (or a keystroke to continue): '
    read,line
    split=strsplit(line,' ',/extract)
    if (n_elements(split) eq 2) then range=float(split)
  end

  'xrange'  : begin
   if (n_elements(xrange) eq 2) then begin
      print,'Display range for wavelength is spectra is ',$
        xrange[0],' to ',xrange[1]
   endif else print,'Display range for wavelength is spectra is ', xrange[0]
    print,'Enter a new (space-delimited) range (or a keystroke to continue): '
    read,line
    split=strsplit(line,' ',/extract)
    if (n_elements(split) eq 2) then xrange=float(split)
  end

  'filter' : begin
    print,'Filter toggled.'
    case filter of
      0 : filter=1
      else : filter=0
    endcase
  end

  'zoom'   : begin
    print,'Zoom factor for image size is ',zoom
    print,'Enter new zoom (or a non-numeric key to continue): '
    zq=byte(get_kbrd(1))
    if (zq ge 49 and zq le 39) then zoom=zq-47
  end

  'help'   : begin
    print," "
    print,"b u k         - Move up one line in file"
    print,"r             - Display and/or modify data range"
    print,"x             - Display and/or modify wavelength range (for spectra)"
    print,"z             - Display and/or modify zoom factor for image size"
    print,"t             - to move to top of file (first line)"
    print,"l             - to move to last line of file (bottom)"
    print,"f             - to change data filtering (if implemented)"
    print,"q e ESC       - To quit"
    print,"h             - For help (this screen)"
    print,"~             - don't hit that key!"
    print,"anything else - Move down one line in file"
  end

  'bad'    : begin
    print,"You were asked not to hit the '~'."
  end

endcase


END
