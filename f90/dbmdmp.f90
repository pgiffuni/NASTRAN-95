SUBROUTINE dbmdmp
!********************************************************************
!     DBMDMP - DUMPS THE IN MEMORY DATA BASE DIRECTORY
!********************************************************************
 INCLUDE 'DSIOF.COM'
 COMMON / zzzzzz / mem(4)
 COMMON / system / isysbf, iwr
 WRITE ( iwr, 900 ) idbbas, idbfre, idbdir, indbas, indclr, indcbp  &
     ,                  nblock, lenalc, iocode, ifilex, NAME,   maxalc  &
     ,                  maxblk, maxdsk, idblen, idbadr, ibasbf, inddir  &
     ,                  numopn, numcls, numwri, numrea, lenopc
 900   FORMAT(/,' CONTENTS OF / DBM / FOLLOW:'  &
     ,/,' IDBBAS =',i8,' IDBFRE =',i8,' IDBDIR =',i8,' INDBAS =',i8  &
     ,/,' INDCLR =',i8,' INDCBP =',i8,' NBLOCK =',i8,' LENALC =',i8  &
     ,/,' IOCODE =',i8,' IFILEX =',i8,' NAME   =',i8,' MAXALC =',i8  &
     ,/,' MAXBLK =',i8,' MAXDSK =',i8,' IDBLEN =',i8,' IDBADR =',i8  &
     ,/,' IBASBF =',i8,' INDDIR =',i8,' NUMOPN =',i8,' NUMCLS =',i8  &
     ,/,' NUMWRI =',i8,' NUMREA -',i8,' LENOPC =',i8 )
 WRITE ( iwr, 901 )
 901   FORMAT(/,' CONTENTS OF FCB FOLLOW:',/)
 DO  i = 1, 80
   WRITE ( iwr, 902 ) i, ( fcb(k,i),k=1,15)
   902   FORMAT(i3,'-',i3,i7,4I5,i12,i2,4I7,2A4,i4)
 END DO
 CALL dbmdia
!      WRITE ( IWR, 906 )
!      WRITE ( IWR, 907 )
 next = idbfre
 itotal = 0
 itotbk = 0
 icnt = 0
 IF ( next == 0 ) GO TO 40
 30    icnt = icnt + 1
 IF ( next == 0 ) GO TO 50
 ival = next
 ivalp= mem(next)
 ivaln= mem(next+1)
 IF ( mem(next  ) == 0 ) ivalp = 0
 IF ( mem(next+1) == 0 ) ivaln = 0
 itotal = itotal + mem(next+2)
 itotbk = itotbk + 1
!      WRITE ( IWR, 908 ) ICNT,IVAL,IVALP,IVALN,MEM(NEXT+2)
 next = mem( next+1 )
 GO TO 30
 40    CONTINUE
!      WRITE( IWR, 909 )
 GO TO 60
 50    CONTINUE
!     WRITE( IWR, 910 ) ITOTAL, ITOTBK
 60    CONTINUE
 700   RETURN
 906   FORMAT(///,31X,' DUMP OF FREE CHAIN',/  &
     ,13X,' ( BLOCK ADDRESSES IN BYTES,  BLOCK LENGTHS IN WORDS )',/)
 907   FORMAT(10X,  &
     '  BLOCK NO    BLOCK ADDRESS  PREV. BLOCK   NEXT BLOCK    LENGTH')
 908   FORMAT( i17,i20,i13,i13,i10)
 909   FORMAT(//' *************** NO FREE SPACE REMAINS **************')
 910   FORMAT(///,' TOTAL FREE SPACE IN WORDS            =',i10  &
     ,/,        ' NUMBER OF BLOCKS IN FREE SPACE CHAIN =',i10)
END SUBROUTINE dbmdmp
