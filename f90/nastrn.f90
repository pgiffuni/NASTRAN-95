PROGRAM nastrn
     
 CHARACTER (LEN=80) :: value
 CHARACTER (LEN=5) :: tmp
 INTEGER :: sperlk
 REAL :: systm(94)
 COMMON / lstadd / lastad
 COMMON / system / isystm(94),sperlk
 COMMON / sofdsn / sdsn(10)
 COMMON / logout / lout
 COMMON / resdic / irdict, iropen
 COMMON / zzzzzz / iz(400000000)
 COMMON / dbm    / idbbas, idbfre, idbdir, indbas, indclr, indcbp  &
     ,                 nblock, lenalc, iocode, ifilex, NAME,   maxalc  &
     ,                 maxblk, maxdsk, idblen, idbadr, ibasbf, inddir  &
     ,                 numopn, numcls, numwri, numrea, lenopc
 INCLUDE 'NASNAMES.COM'
 CHARACTER (LEN=80) :: sdsn
 EQUIVALENCE    ( isystm, systm )
 lenopc = 400000000
 
!     SAVE STARTING CPU TIME AND WALL CLOCK TIME IN /SYSTEM/
 
 isystm(18) = 0
 CALL second (systm(18))
 CALL waltim (isystm(32))
 
!     EXECUTE NASTRAN SUPER LINK
 
 LEN = 80
 value = ' '
 CALL btstrp
 CALL getenv ( 'DBMEM', value )
 READ ( value, * ) idblen
 CALL getenv ( 'OCMEM', value )
 READ ( value, * ) iocmem
 IF ( iocmem <= lenopc ) GO TO 10
 PRINT *,' LARGEST VALUE FOR OPEN CORE ALLOWED IS:',lenopc
 CALL mesage ( -61, 0, 0 )
 10    IF ( idblen /= 0 ) idblen = lenopc - iocmem
 lastad = locfx( iz( iocmem ) )
 IF ( idblen /= 0 ) idbadr = locfx( iz( iocmem+1 ) )
 lenopc = iocmem
 CALL dbmint
 lout   = 3
 irdict = 4
 sperlk = 1
 isystm(11) = 1
 value = ' '
 CALL getenv ( 'RFDIR',  rfdir  )
 value = ' '
 CALL getenv ( 'DIRCTY', dirtry )
 LEN = INDEX( dirtry, ' ' ) - 1
 DO  i = 1, 90
   IF ( i <= 9 ) WRITE ( tmp, 901 ) i
   IF ( i > 9 ) WRITE ( tmp, 902 ) i
   901   FORMAT('scr',i1)
   902   FORMAT('scr',i2)
   dsnames( i ) = dirtry(1:LEN)//'/'//tmp
 END DO
 CALL getenv ( 'LOGNM', LOG )
 dsnames(3) = LOG
 CALL getenv ( 'OPTPNM', optp )
 dsnames(7)  = optp
 CALL getenv ( 'NPTPNM', nptp )
 dsnames(8)  = nptp
 CALL getenv ( 'FTN11', out11 )
 dsnames(11) = out11
 CALL getenv ( 'FTN12', in12 )
 dsnames(12) = in12
 CALL getenv ( 'FTN13', value )
 dsnames(13) = value
 CALL getenv ( 'FTN14', value )
 dsnames(14) = value
 CALL getenv ( 'FTN15', value )
 dsnames(15) = value
 CALL getenv ( 'FTN16', value )
 dsnames(16) = value
 CALL getenv ( 'FTN17', value )
 dsnames(17) = value
 CALL getenv ( 'FTN18', value )
 dsnames(18) = value
 CALL getenv ( 'FTN19', value )
 dsnames(19) = value
 CALL getenv ( 'FTN20', value )
 dsnames(20) = value
 CALL getenv ( 'FTN21', value )
 dsnames(21) = value
 CALL getenv ( 'PLTNM', plot )
 dsnames(10) = plot
 CALL getenv ( 'DICTNM', dic )
 dsnames(4) = dic
 CALL getenv ( 'PUNCHNM', punch )
 dsnames(1) = punch
 CALL getenv ( 'SOF1', value )
 sdsn(1) = value
 CALL getenv ( 'SOF2', value )
 sdsn(2) = value
 CALL getenv ( 'SOF3', value )
 sdsn(3) = value
 CALL getenv ( 'SOF4', value )
 sdsn(4) = value
 CALL getenv ( 'SOF5', value )
 sdsn(5) = value
 CALL getenv ( 'SOF6', value )
 sdsn(6) = value
 CALL getenv ( 'SOF7', value )
 sdsn(7) = value
 CALL getenv ( 'SOF8', value )
 sdsn(8) = value
 CALL getenv ( 'SOF9', value )
 sdsn(9) = value
 CALL getenv ( 'SOF10', value )
 sdsn(10) = value
 OPEN (  3, FILE=dsnames(3) ,STATUS='UNKNOWN')
 IF ( dsnames(11) /= 'none' ) OPEN ( 11, FILE=dsnames(11),STATUS='UNKNOWN')
 IF ( dsnames(12) /= 'none' ) OPEN ( 12, FILE=dsnames(12),STATUS='UNKNOWN')
 IF ( dsnames(10) /= 'none' ) OPEN ( 10, FILE=dsnames(10),STATUS='UNKNOWN')
 IF ( dsnames(4) /= 'none' ) OPEN ( 4, FILE=dsnames(4),STATUS='UNKNOWN')
 IF ( dsnames(1) /= 'none' ) OPEN ( 1, FILE=dsnames(1),STATUS='UNKNOWN')
 CALL xsem00
 STOP
END PROGRAM nastrn
