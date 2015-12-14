SUBROUTINE dbmfdp
!********************************************************************
!     DBMFDP- DUMPS THE DIRECTORY CHAIN OF A GIVEN FILE.
!             ARGUMENT IDIR IS THE IN-MEMORY DIRECTORY FOR THE FILE
!********************************************************************
 COMMON / system / isysbf, iwr
 INCLUDE  'ZZZZZZ.COM'
 INCLUDE  'DSIOF.COM'
 ibase  = locfx( mem )
 ival2  = ibase + fcb(  9, ifilex )
 ival3  = ibase + fcb( 10, ifilex )
 ival4  = ibase + fcb( 11, ifilex )
 INDEX  = fcb( 10, ifilex )
 lblock = mem( INDEX+3 )
 WRITE ( iwr, 902 ) ifilex, ival2, ival3, ival4, fcb(12,ifilex)
 WRITE ( iwr, 903 )
 next = fcb(  9, ifilex )
 icnt = 0
 IF ( next == 0 ) GO TO 25
 20    icnt = icnt + 1
 IF ( next == 0 ) GO TO 30
 ival = ibase + next
 ivalp= ibase + mem(next)
 ivaln= ibase + mem(next+1)
 IF ( mem( next   ) == 0 ) ivalp = 0
 IF ( mem( next+1 ) == 0 ) ivaln = 0
 WRITE ( iwr, 904 ) mem(next+3),mem(next+7),ival,ivalp,ivaln,mem(next+2)
 990   FORMAT( 12(8(1X,i8),/))
 next = mem( next+1 )
 GO TO 20
 25    WRITE( iwr, 907 )
 30    CONTINUE
 WRITE( iwr, 908 )
 RETURN
 902   FORMAT(///,25X,' DUMP OF FILE CHAIN FOR UNIT=',i6,/  &
     ,14X,'( BLOCK ADDRESSES ARE IN WORDS,  BLOCK LENGTHS IN WORDS)',/ ,/,7X,  &
     ' FIRST BLOCK ADDRESS   ',i12,'   LAST BLOCK ADDRESS      ',i12 ,/,7X,  &
     ' CURRENT BLOCK ADDRESS ',i12,'   ORIGINAL BUFFER ADDRESS ',i12)
 903   FORMAT(/, '  IN-MEM     BUFFER',/  &
     ,' BLOCK NO.  BLOCK NO  BLOCK ADDRESS  PREV. BLOCK   NEXT BLOCK ' ,' LENGTH')
 904   FORMAT( i9,i11,5X,i12,7X,i12,5X,i12,i12)
 907   FORMAT(//' *************** NO BLOCK ALLOCATED TO FILE **********')
 908   FORMAT(///)
END SUBROUTINE dbmfdp
