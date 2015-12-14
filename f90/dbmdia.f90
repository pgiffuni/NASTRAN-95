SUBROUTINE dbmdia
!********************************************************************
!     DBMDIA - DUMPS THE IN MEMORY DATA BASE DIRECTORY
!********************************************************************
 INCLUDE  'DSIOF.COM'
 INCLUDE  'ZZZZZZ.COM'
 COMMON / system / isysbf, iwr
 INTEGER :: SCRATCH(2)
 DATA              SCRATCH / 4HSCRA , 4HTCHX /
 iblksz = isysbf - 4
 itoti  = 0
 itotx  = 0
 WRITE ( iwr, 903 )
 DO  i = 1, 80
   IF ( i == 7 ) CYCLE
   IF ( fcb( 9,i ) == 0 .AND. fcb( 5,i ) == 0 ) CYCLE
   INDEX = fcb( 10, i )
   iintb = 0
   iextb = 0
   IF ( fcb( 9,i ) /= 0 ) iintb = mem( INDEX+3 )
   itoti = itoti + iintb
   IF ( fcb( 5,i ) /= 0 ) iextb = fcb(6,i) - fcb( 5,i) + 1
   IF ( iextb >= fcb( 7, ifilex ) ) CYCLE
   itotx = itotx + iextb
   IF ( fcb( 13,i ) /= 0 ) GO TO 15
   fcb( 13,i ) = SCRATCH(1)
   fcb( 14,i ) = SCRATCH(2)
   15    CONTINUE
   WRITE ( iwr, 904 ) i, fcb( 13,i ), fcb( 14,i ), fcb( 4,i )  &
       ,                  iintb, iextb
 END DO
 WRITE ( iwr, 905 ) itoti, itotx
!      WRITE ( IWR, 906 ) MAXBLK, MAXDSK, MAXALC, IBLKSZ
 700   RETURN
 903   FORMAT(///,27X,' MEMORY DATA BASE DIRECTORY',//,  &
     '    UNIT    NAME   CURRENT  IN-MEM' ,'   DISK ',/,  &
     '                    BLOCK   BLOCKS' ,'  BLOCKS ',/)
 904   FORMAT(i7,3X,2A4,2X,i6,2X,i6,2X,i6 )
 905   FORMAT(/,' CURRENT IN-MEMORY BLOCKS =',i16  &
     ,/,' CURRENT DISK BLOCKS      =',i16 )
 906   FORMAT(/,' MAXIMUM IN-MEMORY BLOCKS USED                   =',i16  &
     ,/,' MAXIMUM DISK BLOCKS WRITTEN                     =',i16  &
     ,/,' BLOCKS INITIALLY ALLOCATED FOR THE IN-MEMORY DB =',i16  &
     ,/,' BLOCK SIZE                                      =',i16)
END SUBROUTINE dbmdia
