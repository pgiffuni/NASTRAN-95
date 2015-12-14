SUBROUTINE dbmdfc
!********************************************************************
!     DBMDFC - DUMPS THE FREE CHAIN
!********************************************************************
 INCLUDE  'DSIOF.COM'
 INCLUDE  'ZZZZZZ.COM'
 COMMON / system / isysbf, iwr
 WRITE ( iwr, 906 )
 WRITE ( iwr, 907 )
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
 WRITE ( iwr, 908 ) icnt,ivalp,ival,ivaln,mem(next+2)
 next = mem( next+1 )
 GO TO 30
 40    CONTINUE
 WRITE( iwr, 909 )
 GO TO 60
 50    CONTINUE
 WRITE( iwr, 910 ) itotal, itotbk
 60    CONTINUE
 700   RETURN
 906   FORMAT(///,31X,' DUMP OF FREE CHAIN',/  &
     ,13X,' ( BLOCK ADDRESSES IN WORDS,  BLOCK LENGTHS IN WORDS )',/)
 907   FORMAT(10X,  &
     '  BLOCK NO    PREV. BLOCK    BLOCK ADDRESS NEXT BLOCK    LENGTH')
 908   FORMAT( i17,i20,i13,i13,i10)
 909   FORMAT(//' *************** NO FREE SPACE REMAINS **************')
 910   FORMAT(///,' TOTAL FREE SPACE IN WORDS            =',i10  &
     ,/,        ' NUMBER OF BLOCKS IN FREE SPACE CHAIN =',i10)
END SUBROUTINE dbmdfc
