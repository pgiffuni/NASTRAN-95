SUBROUTINE dbmmgr ( opcode )
!*********************************************************************
!        / FCB /
!            FCB(1,I) - OPEN FLAG
!            FCB(2,I) - BUFFER ADDRESS
!            FCB(3,I) - CURRENT CLR
!            FCB(4,I) - CURRENT BLOCK NUMBER
!            FCB(5,I) - FIRST BLOCK NUMBER WRITTEN TO THIS FILE
!            FCB(6,I) - LAST BLOCK NUMBER WRITTEN TO THIS FILE
!            FCB(7,I) - MAXIMUM NUMBER OF BLOCKS TO BE ALLOCATED
!                       TO THIS FILE
!            FCB(8,I) - =0, IF NO MATRIX STRINGS WRITTEN TO FILE
!                       =1, OTHERWISE, USED TO INITIALIZE COLUMN
!                           NUMBER TO 1.
!            FCB(9,I) - INDEX TO FIRST IN-MEMORY BLOCK
!            FCB(10,I)- INDEX TO LAST IN-MEMORY BLOCK
!            FCB(11,I)- INDEX TO CURRENT IN-MEMORY BLOCK
!            FCB(12,I)- ORIGINAL BUFFER ADDRESS
!            FCB(13-14,I) - DMAP FILE NAME (2A4)
!            FCB(15,I)- OPEN FLAG FOR EXTERNAL FILE
!        / DBM/
!            IDBBAS - (INPUT)-INDEX TO IN-MEMORY DATA BASE RELATIVE
!                              TO /DBM/
!            IDBFRE - (INPUT)-INDEX TO FREE CHAIN OF IN-MEMORY DATA
!                              BASE RELATIVE TO /DBM/
!            IDBDIR - (INPUT)-INDEX TO FIRST DIRECTORY BLOCK
!            MAXALC - (OUTPUT)-MAXIMUM NUMBER OF BLOCKS AVAILABLE FOR
!                              JOB
!            MAXBLK - (OUTPUT)-MAXIMUM NUMBER OF BLOCKS ALLOCATED(JOB)
!            MAXDSK - (OUTPUT)-MAXIMUM NUMBER OF BLOCKS WRITTEN TO
!                              TO DISK
!            LENALC - (OUTPUT)-LENGTH OF EACH ALLOCATED BLOCK
!            IOCODE - (INPUT) -IO-CODE FOR OPEN/CLOSE CALL
!            IFILEX - (INPUT) -FILE NUMBER FOR GINO FILE IN /XFIAT/
!            NBLOCK - (INPUT/OUTPUT) -BLOCK NUMBER BEING REFERENCED
!            NAME   - (INPUT) -GINO FILE NAME (E.G., 101,201,303,...)
!            INDBAS - INDEX TO START OF BUFFER RELATIVE TO /ZZZZZZ/
!            INDCLR - INDEX TO CLR WITHIN BUFFER RELATIVE TO /ZZZZZZ/
!            INDCBP - INDEX TO CBP WITHIN BUFFER RELATIVE TO /ZZZZZZ/
!        FREE CHAIN FORMAT (ALSO, ALL BLOCKS ALLOCATED)
!               IDBFRE==> WORD 0    POINTER TO PREVIOUS FREE BLOCK
!                                      IN CHAIN, ALWAYS 0 FOR 1ST BLK)
!                         WORD 1    POINTER TO NEXT BLOCK IN CHAIN
!                                      -INITIALLY SET TO ZERO)
!                         WORD 2    NUMBER OF FREE WORDS IN BLOCK
!                         WORD 3    RELATIVE BLOCK NUMBER
 
!         OPCODE
!           1    OPEN
!                  /GINOX/ IOCODE = 0 ; READ WITH REWIND
!                                 = 1 ; WRITE WITH REWIND
!                                 = 2 ; READ WITHOUT REWIND
!                                 = 3 ; WRITE WITHOUT REWIND
!           2    CLOSE
!                  /GINOX/ IOCODE = 1 ; CLOSE WITH REWIND
!                                    (OTHERWISE NO REWIND)
!           3    REWIND
!           4    WRITE
!           5    READ
!           6    POSITION FILE
!                  NBLOCK = BLOCK NUMBER TO POSITION TO
!           7    DELETE FILE
!           8    PROCESS WRTBLK REQUEST (SUBSTRUCTURING)
!           9    PROCESS RDBLK REQUEST (SUBSTRUCTURING)
!********************************************************************
 
 INTEGER, INTENT(IN)                      :: opcode
 
 INTEGER :: case  / 4HCASE /
 INTEGER :: xycd  / 4HXYCD /
 INTEGER :: pcdb  / 4HPCDB /
 INTEGER :: pool  / 4HPOOL /
 INTEGER :: xpdt  / 4HXPDT /
 INCLUDE    'DSIOF.COM'
 COMMON / xfist  / fist(10)
 COMMON / xfiat  / fiat(10)
 COMMON / zzzzzz / mem(4)
 COMMON / system / isysbf, iwr
 DATA     lenbuf / 0 /
 
 IF ( lenbuf /= 0 ) GO TO 10
! SET UP BLOCK ALLOCATIONS FOR DOUBLE WORD BOUNDARIES
 ibasbf = locfx( mem )
 lenbuf = isysbf - 3 + 8
 lenalc = lenbuf
 nbuff3 = isysbf - 4
 itest  = MOD( lenbuf,2)
 IF ( itest /= 0 ) lenbuf = lenbuf + 1
 10    IF ( idbdir /= 0 ) GO TO 30
! OPCODES OF 8 AND 9 HAVE NO PURPOSE WHEN THERE IS NO USE OF THE
! IN-MEMORY DATA BASE
 IF ( opcode == 8 .OR. opcode == 9 ) GO TO 7777
!  CALL DBMIO DIRECTLY, NO IN-MEMORY DATA BASE
 20    CALL dbmio ( opcode )
 GO TO 7777
 30    IF ( NAME > 100 .AND. NAME < 400 ) GO TO 50
!30    IF ( NAME .GT. 300 .AND. NAME .LT. 400 ) GO TO 50
!  CHECK FOR CASECC, XYCD, AND PCDB (SETUP IN FIAT FOR PREFACE)
 IF ( NAME == case ) GO TO 50
 IF ( NAME == xycd ) GO TO 50
 IF ( NAME == pcdb ) GO TO 50
 IF ( NAME == xpdt ) GO TO 50
 IF ( NAME == pool ) GO TO 50
! OPCODES OF 8 AND 9 HAVE NO PURPOSE WHEN THERE IS NO USE OF THE
! IN-MEMORY DATA BASE
 IF ( opcode == 8 .OR. opcode == 9 ) GO TO 7777
! CALL DBMIO DIRECTLY BECAUSE THIS IS AN EXECUTIVE FILE
 IF ( fcb(  9, ifilex ) /= 0 ) CALL dbmrel
 GO TO 20
 50    CONTINUE
!      IF ( IFILEX .NE. 48 ) GO TO 55
!      IF ( NAME .NE. 307 ) GO TO 55
!      WRITE(IWR,40646)OPCODE,IOCODE,NBLOCK,IFILEX,NAME,INDBAS
 40646 FORMAT(/,' OPCODE,IOCODE,NBLOCK,IFILEX,NAME,INDBAS=',6I6)
!      WRITE(IWR,40647)(MEM(INDBAS+KB),KB=-4,20)
 40647 FORMAT(' INPUT BUFFER HAS=',/,10(4(1X,z8),/))
!      WRITE(6,44770)(FCB(K,IFILEX),K=1,15)
 44770 FORMAT(' ENTERRED FCB=',/,2(5I8,/),2I8,4X,2A4,4X,i8)
!      CALL DBMFDP
 55    CONTINUE
 SELECT CASE ( opcode )
   CASE (    1)
     GO TO  100
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 300
   CASE (    4)
     GO TO 400
   CASE (    5)
     GO TO 500
   CASE (    6)
     GO TO 600
   CASE (    7)
     GO TO 700
   CASE (    8)
     GO TO 800
   CASE (    9)
     GO TO 900
 END SELECT
!****************
! OPEN CODE *********************************************************
!****************
 100   CONTINUE
 fcb(  1, ifilex ) = iocode
 fcb( 12, ifilex ) = fcb(  2, ifilex )
 IF ( fcb( 9, ifilex ) /= 0 ) GO TO 130
! CHECK TO SEE IF FILE IS SELF CONTAINED ON DISK
 IF ( fcb( 5, ifilex ) /= 0 ) GO TO 120
 105   CONTINUE
 IF ( iocode /= 0 .AND. iocode /= 2 ) GO TO 108
 WRITE ( iwr, 9900 ) ifilex, fcb( 13, ifilex), fcb( 14, ifilex )
 9900  FORMAT(///,' DBMMGR ERROR, ATTEMPT TO OPEN FOR READ OR WRITE APP'  &
     ,'END:' ,/,' UNIT-',i4,'  NAME=',2A4,' WHICH DOES NOT EXIST.')
!      CALL DBMDMP
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
 108   CONTINUE
! NEW FILE NAME FOR IFILEX, RELEASE ANY PREVIOUSLY ALLOCATED BLOCKS
 IF ( fcb(  9, ifilex ) /= 0 ) CALL dbmrel
! CREATE FILE ENTRY IN FCB
 DO  i = 3,11
   IF ( i == 7 ) CYCLE
   fcb(  i, ifilex ) = 0
 END DO
 fcb(  4, ifilex ) = 1
 nblock            = 1
 115   CONTINUE
! ALLOCATE FIRST BLOCK
 CALL dbmalb ( lenbuf, nexblk )
 IF ( nexblk <= 0 ) GO TO 120
 fcb(  9, ifilex ) = nexblk
 fcb( 10, ifilex ) = nexblk
 fcb( 11, ifilex ) = nexblk
! INITIALIZE PREVIOUS, NEXT, LENGTH AND BLOCK NUMBER FOR ALLOCATED BLK
 mem( nexblk   ) = 0
 mem( nexblk+1 ) = 0
 mem( nexblk+2 ) = lenbuf
 mem( nexblk+3 ) = 1
 fcb(  2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 4)
 GO TO 7000
! NO MORE SPACE WITHIN IN-MEMORY DATA BASE, USE I/O
 120   CALL dbmio ( opcode )
 GO TO 7777
! FILE EXISTS IN IN-MEMORY DATA BASE
 130   CONTINUE
 IF ( iocode == 0 ) GO TO 150
 IF ( iocode == 1 ) GO TO 160
 IF ( iocode == 2 ) GO TO 170
 IF ( iocode == 3 ) GO TO 180
! FILE IS OPENED FOR READ WITH REWIND
 150   CONTINUE
 nexblk = fcb( 9, ifilex )
 IF ( nexblk > 0 ) GO TO 155
 WRITE ( iwr, 9910 ) ifilex
 9910  FORMAT(///,' DBMMGR ERROR, ATTEMPT TO READ FILE WITH NO BLOCKS'  &
     /,' UNIT=',i4)
!      CALL DBMDMP
 CALL dsmsg ( 777 )
 CALL mesage( -61, 0, 0 )
 155   CONTINUE
 fcb( 11, ifilex ) = nexblk
 fcb(  4, ifilex ) = 1
 nblock            = 1
 fcb(  2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 3 )
 GO TO 7000
! FILE IS OPENED FOR WRITE WITH REWIND
 160   CONTINUE
 GO TO 105
! FILE IS OPENED FOR READ WITHOUT REWIND
 170   CONTINUE
 nexblk = fcb(  10, ifilex )
 lastib = mem( nexblk+3 )
 nblock = fcb( 4, ifilex )
 IF ( fcb( 4, ifilex ) > lastib ) GO TO 120
 IF ( fcb( 4, ifilex ) == 1 ) GO TO 150
 nexblk = fcb( 11, ifilex )
 iblk1 = fcb(  4, ifilex )
 iblk2 = mem( nexblk+3 )
 iblk3 = mem( nexblk+7 )
 fcb(  2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
! CHECK THAT CURRENT BLOCK NUMBER MATCHES BLOCK NO. IN IN-MEM BLK
 IF ( iblk1 == iblk2 .AND. iblk1 == iblk3 ) GO TO 7000
 GO TO 190
! FILE IS OPENED FOR WRITE WITHOUT REWIND
 180   CONTINUE
 nexblk = fcb(  10, ifilex )
 lastib = mem( nexblk+3 )
 IF ( fcb( 4, ifilex ) > lastib ) GO TO 120
!======      IF ( FCB( 4, IFILEX ) .EQ. 1      ) GO TO 160
 nexblk = fcb( 11, ifilex )
! IGNORE ANY PREVIOUSLY WRITTEN BLOCKS FOR THIS FILE
 fcb(  5, ifilex ) = 0
 fcb(  6, ifilex ) = 0
 iblk1  = fcb(  4, ifilex )
 iblk2  = mem( nexblk+3 )
 iblk3  = mem( nexblk+7 )
 fcb(  2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
! CHECK THAT CURRENT BLOCK NUMBER MATCHES BLOCK NO. IN IN-MEM BLK
 IF ( iblk1 == iblk2 .AND. iblk1 == iblk3 ) GO TO 7000
 190   CONTINUE
 WRITE ( iwr, 9911 ) ifilex, iblk1, iblk2, iblk3
 9911  FORMAT(///' BLOCK NUMBERS INCONSISTANT ON OPEN IN DBMMGR'  &
     ,/,' UNIT =',i4 ,/,' BLOCK NUMBER EXPECTED (IN FCB)  =',i8  &
     ,/,' BLOCK NUMBER IN IN-MEMORY BLOCK =',i8  &
     ,/,' BLOCK NUMBER IN BUFFER          =',i8 )
!      CALL DBMDMP
 CALL dbmfdp
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
!****************
! CLOSE CODE ********************************************************
!****************
 200   CONTINUE
! CHECK TO SEE IF FILE HAS IN-MEMORY BLOCKS
 IF ( fcb(  9, ifilex ) /= 0 ) GO TO 220
 210   CALL dbmio ( opcode )
 GO TO 7000
 220   CONTINUE
!WKBDB SPR94012 10/94
!      IF ( IOCODE .NE. 1 ) GO TO 225
!C CLOSE FILE WITH REWIND
!      FCB( 11, IFILEX ) = FCB(  9, IFILEX )
!      FCB(  4, IFILEX ) = 1
!      IF ( FCB( 5, IFILEX ) .NE. 0 ) GO TO 210
!WKBDE SPR94012 10/94
! IF FILE IS OPENED FOR READ THAN GO COMPUTE STATISTICS
 225   IF ( fcb(  1, ifilex ) == 0.OR. fcb(  1, ifilex ) == 2 ) GO TO 240
 IF ( fcb( 15, ifilex ) /= 0 ) GO TO 240
! FILE OPENED FOR WRITE AND FILE NOT SPILLED TO DISK, THEN
! RELEASE LAST ALLOCATED BLOCK, BECAUSE IT WAS NOT USED
 nexblk = fcb( 11, ifilex )
! RESET LAST BLOCK POINTER, GET PREVIOUS BLOCK ALLOCATED
!WKBNB SPR94012 10/94
 228   iblock = mem( nexblk+3 )
! CHECK IF LAST BLOCK NOT USED, THERE COULD HAVE BEEN A BACKPSPACE BACK
! TO A PREVIOUS USED BLOCK (CAUSED BY CLOSE CALLING DSBRC1 TO BACKSPACE
! OVER AN EOF THAT WAS AT THE END OF A PREVIOUS BLOCK).
 IF ( iblock > nblock ) GO TO 230
 nexblk = mem( nexblk+1 )
 IF ( nexblk == 0 ) GO TO 240
 GO TO 228
 230   CONTINUE
!WKBNE SPR94012 10/94
 indblk            = mem( nexblk )
 fcb( 10, ifilex ) = indblk
 fcb( 11, ifilex ) = indblk
 fcb(  4, ifilex ) = mem( indblk+3 )
 fcb(  2, ifilex ) = locfx( mem( indblk+4 ) ) - ibasbf + 1
 CALL dbmrlb( nexblk )
!WKBNB SPR94012 10/94
 240   IF ( iocode /= 1 ) GO TO 245
! CLOSE FILE WITH REWIND
 fcb( 11, ifilex ) = fcb(  9, ifilex )
 fcb(  4, ifilex ) = 1
!WKBNE SPR94012 10/94
!WKBR  SPR94012 10/94
!240   IF ( FCB( 5, IFILEX ) .NE. 0 ) CALL DBMIO ( OPCODE )
 245   IF ( fcb( 5, ifilex ) /= 0 ) CALL dbmio ( opcode )
 IF ( fcb( 5, ifilex ) <= fcb( 6, ifilex ) ) GO TO 7000
! SPECIAL CASE, LAST BLOCK ALLOCATED WAS FOR DISK BUT NEVER USED, RESET
! INDBAS BACK TO LAST IN-MEMORY BLOCK
 nexblk = fcb( 10, ifilex )
 fcb( 2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 fcb( 5, ifilex ) = 0
 fcb( 6, ifilex ) = 0
 fcb(11, ifilex ) = fcb( 10, ifilex )
 GO TO 7000
!****************
! REWIND OPCODE *****************************************************
!****************
 300   CONTINUE
! IF FILE IS ON EXTERNAL FILE CALL DBMIO DIRECTLY
 IF ( fcb( 9, ifilex ) /= 0 ) GO TO 320
 CALL dbmio ( opcode )
 GO TO 7777
 320   CONTINUE
 nexblk            = fcb(  9, ifilex )
 fcb( 11, ifilex ) = nexblk
 fcb(  4, ifilex ) = 1
! REPLACE BUFFER ADDRESS IN FCB
 fcb( 2,ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 3 )
 iocode = 0
 IF ( fcb( 5, ifilex ) /= 0 ) CALL dbmio ( 2 )
 GO TO 7000
!****************
! WRITE CODE ********************************************************
!****************
 400   CONTINUE
! CHECK TO SEE IF THIS BLOCK IS ON EXTERNAL FILE
 IF ( fcb( 15, ifilex ) /= 0 ) GO TO 450
! CHECK THAT BLOCK NUMBER MATCHES
 nexblk = fcb( 11, ifilex )
 iblk1  = fcb(  4, ifilex )
 iblk2  = mem( nexblk+3 )
 iblk3  = mem( nexblk+7 )
 IF ( iblk1 == iblk2 .AND. iblk1 == iblk3 ) GO TO 410
 WRITE ( iwr, 9940 ) ifilex, iblk1, iblk2, iblk3
 9940  FORMAT(///' BLOCK NUMBERS INCONSISTANT ON WRITE IN DBMMGR'  &
     ,/,' UNIT = ',i4 ,/,' BLOCK NUMBER EXPECTED (IN FCB)  =',i8  &
     ,/,' BLOCK NUMBER IN IN-MEMORY BLOCK =',i8  &
     ,/,' BLOCK NUMBER IN BUFFER          =',i8 )
!      CALL DBMDMP
 CALL dbmfdp
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
 410   CONTINUE
 fcb(  4, ifilex ) = fcb(  4, ifilex ) + 1
 nexblk = mem( indbas-3 )
 IF ( nexblk == 0 ) GO TO 420
! USE EXISTING BLOCK ALREADY ALLOCATED FROM PREVIOUS OPEN FOR WRITE
 fcb( 11, ifilex) = nexblk
 fcb( 2,ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 4 )
 GO TO 7000
 420   CONTINUE
 CALL dbmalb ( lenbuf, nexblk )
 IF ( nexblk <= 0 ) GO TO 440
! ANOTHER BLOCK SUCCESSFULLY ALLOCATED, CONNECT TO CHAIN
 indblk          = fcb( 11, ifilex )
 mem( indblk+1 ) = nexblk
 mem( nexblk   ) = indblk
 mem( nexblk+1 ) = 0
 mem( nexblk+2 ) = lenbuf
 mem( nexblk+3 ) = fcb(  4, ifilex )
 fcb( 10, ifilex) = nexblk
 fcb( 11, ifilex) = nexblk
 fcb( 2,ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 4 )
 GO TO 7000
! NO MORE SPACE IN IN-MEMORY DATA BASE, WRITE DATA TO FILE
 440   CONTINUE
! CALL DBMIO TO OPEN EXTERNAL FILE WITH REWIND
 isave  = iocode
 isaveb = nblock
 iocode = 1
 nblock = fcb( 4, ifilex )
 iprblk = indbas
! RESET BUFFER ADDRESS TO BUFFER IN USER'S OPEN CORE
 fcb( 2,ifilex ) = fcb( 12, ifilex )
 indbas = fcb( 2, ifilex )
 CALL dbmio ( 1 )
 iocode   = isave
 nblock   = isaveb
!      WRITE(6,88771)(MEM(IPRBLK+K),K=-4,4)
 88771 FORMAT(' MEMPRBLK=',9(1X,z8))
!      WRITE(6,88772)(MEM(INDBAS+K),K=-4,4)
 88772 FORMAT(' MEMINDBAS=',9(1X,z8))
!      PRINT *,' IFILEX,NBLOCK,IPRBLK,INDBAS=',IFILEX,NBLOCK,
!     & IPRBLK,INDBAS
 CALL dbmmov ( iprblk, indbas, 4 )
!      PRINT *,' MEM(IPRBLK=',MEM(IPRBLK)
!      WRITE(6,88771)(MEM(IPRBLK+K),K=-4,4)
!      WRITE(6,88772)(MEM(INDBAS+K),K=-4,4)
 GO TO 7000
 450   CONTINUE
 CALL dbmio ( opcode )
 GO TO 7777
!****************
! READ CODE *********************************************************
!****************
 500   CONTINUE
 IF ( fcb( 5, ifilex ) == 0 ) GO TO 505
 IF ( fcb( 4, ifilex ) >= ( fcb( 5, ifilex ) - 1 ) ) GO TO 540
 505   fcb( 4, ifilex ) = fcb( 4, ifilex ) + 1
 nexblk = mem( indbas-3 )
 IF ( nexblk > 0 ) GO TO 510
 WRITE ( iwr, 9950 ) fcb( 4, ifilex ), ifilex
 9950  FORMAT(///,' ERROR IN DBMMGR DURING READ',/,' EXPECTED ANOTHER '  &
     ,' IN-MEMORY BLOCK FOR BLOCK=',i8,' UNIT=',i3)
!      CALL DBMDMP
 CALL dbmfdp
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
 510   fcb(  2, ifilex ) = locfx( mem( nexblk+4 ) ) - ibasbf + 1
 fcb( 11, ifilex ) = nexblk
 CALL dbmmov ( indbas, nexblk+4, 3 )
 iblk1 = fcb( 4, ifilex)
 iblk2 = mem( nexblk+3 )
 iblk3 = mem( nexblk+7 )
 IF ( iblk1 == iblk2 .AND. iblk1 == iblk3 ) GO TO 7000
 WRITE ( iwr, 9951 ) ifilex, iblk1, iblk2, iblk3
 9951  FORMAT(///' BLOCK NUMBERS INCONSISTANT ON READ IN DBMMGR'  &
     ,/,' UNIT =',i4 ,/,' BLOCK NUMBER  (IN FCB)          =',i8  &
     ,/,' BLOCK NUMBER IN IN-MEMORY BLOCK =',i8  &
     ,/,' BLOCK NUMBER IN BUFFER          =',i8 )
!      CALL DBMDMP
 CALL dbmfdp
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
! BLOCK IS NOT IN MEMORY, CALL DBMIO
 540   CONTINUE
 IF ( fcb( 15, ifilex ) /= 0 ) GO TO 550
 isave  = iocode
 isaveb = nblock
 iocode = 0
 nblock = fcb( 4, ifilex ) + 1
 iprblk = indbas
 indbas = fcb( 12, ifilex )
 fcb( 2, ifilex ) = indbas
 CALL dbmio ( 1 )
 iocode = isave
 nblock = isaveb
 CALL dbmmov ( iprblk, indbas, 3 )
 GO TO 7777
 550   CONTINUE
 IF ( fcb( 4, ifilex ) > fcb( 6, ifilex ) ) GO TO 570
 indbas = fcb( 12, ifilex )
 fcb( 2, ifilex ) = indbas
 CALL dbmio ( opcode )
 GO TO 7777
 570   CONTINUE
 WRITE ( iwr, 9052 ) ifilex
 9052  FORMAT(///,' DBMMGR ERROR, ATTEMPT TO READ BEYOND EOF' ,/' UNIT=',i5)
!      CALL DBMDMP
 CALL dbmfdp
 CALL dsmsg ( 777 )
 CALL mesage ( -61, 0, 0 )
!****************
! POSITION CODE *****************************************************
!****************
 600   CONTINUE
 IF ( fcb( 5, ifilex ) == 0 ) GO TO 605
 IF ( nblock >= fcb( 5, ifilex ) ) GO TO 690
 605   CONTINUE
! BLOCK IS IN THE IN-MEMORY DATA BASE, WALK CHAIN TO CORRECT BLOCK
 ioff    = 1
 nblk    = nblock - 1
 nexblk  = fcb( 9, ifilex )
 IF ( nblock == 1 ) GO TO 670
 icndex = fcb( 11, ifilex )
 IF ( icndex == 0 ) GO TO 610
 nexblk  = icndex
 icblk   = mem( icndex+3 )
 IF ( icblk == nblock ) GO TO 670
 idiff   = nblock - icblk
 nblk    = IABS( idiff )
 IF ( idiff < 0 ) ioff = 0
 610   CONTINUE
 DO  i = 1, nblk
   nexblk  = mem( nexblk+ioff )
 END DO
! SET DIRECTORY ENTRIES FOR THE POSITIONED BLOCK
 670   fcb( 11, ifilex ) = nexblk
 fcb(  4, ifilex ) = nblock
 fcb(  2, ifilex ) = locfx( mem(nexblk+4) ) - ibasbf + 1
 CALL dbmmov ( indbas, nexblk+4, 3 )
 GO TO 7000
 690   CONTINUE
 IF ( fcb( 15, ifilex ) /= 0 ) GO TO 695
 isave  = iocode
 iocode = 0
 iprblk = indbas
 indbas = fcb( 12, ifilex )
 fcb( 2, ifilex ) = indbas
 fcb( 4, ifilex ) = nblock
 CALL dbmio( 1 )
 iocode = isave
 CALL dbmmov ( iprblk, indbas, 3 )
 GO TO 7777
 695   CONTINUE
 fcb( 4, ifilex ) = nblock
 indbas = fcb( 12, ifilex )
 fcb( 2, ifilex ) = indbas
 CALL dbmio ( opcode )
 GO TO 7777
!****************
! DELETE CODE *******************************************************
!****************
 700   CONTINUE
 IF ( fcb( 9, ifilex ) == 0 ) GO TO 710
 CALL dbmrel
 710   CONTINUE
 CALL dbmio ( 7 )
 DO  k = 1,15
   IF ( k == 7 ) CYCLE
   fcb( k, ifilex ) = 0
 END DO
 GO TO 7777
!****************
! WRTBLK CODE *******************************************************
!****************
! SPECIAL ENTRY FOR SUBSTRUCTURING, MOVE DATA FROM OPENCORE BUFFER
! CALLED BY WRTBLK OF GINO
 800   CONTINUE
 IF ( fcb( 15, ifilex ) == 0 ) GO TO 810
! ORIGINAL BUFFER IS BEING USED BY GINO, JUST RETURN
 GO TO 7777
 810   ind1 = fcb(  2, ifilex )
 ind2 = fcb( 12, ifilex )
 ind1 = ind1 + 2
 ind2 = ind2 + 2
!      PRINT *,' DBMMGR,WRTBLK,IND1,IND2,NBUFF3=',IND1,IND2,NBUFF3
!      PRINT *,' DBMMGR,WRTBLK,INDBAS=',INDBAS
!      WRITE(6,44771)(FCB(K,IFILEX),K=1,15)
!      WRITE(6,44772)(MEM(IND2+K),K=1,8)
 44772 FORMAT(' DBMMGR,BUFFER,IND2=',8(1X,z8))
 DO  i = 1, nbuff3
   mem( ind1+i ) = mem( ind2+i )
 END DO
 GO TO 7000
!****************
! RDBLK  CODE *******************************************************
!****************
! SPECIAL ENTRY FOR SUBSTRUCTURING, MOVE DATA TO ORIGINAL BUFFER IF
! THE IN-MEMORY DATA BASE IS BEING USED
! CALLED BY RDBLK
 900   CONTINUE
 IF ( fcb( 15, ifilex ) == 0 ) GO TO 910
! ORIGINAL BUFFER IS BEING USED, JUST RETURN
 GO TO 7777
 910   ind1 = fcb(  2, ifilex )
 ind2 = fcb( 12, ifilex )
 ind1 = ind1 + 2
 ind2 = ind2 + 2
!      PRINT *,' DBMMGR,RDBLK,IND1,IND2,NBUFF3=',IND1,IND2,NBUFF3
!      PRINT *,' DBMMGR,RDBLK,INDBAS=',INDBAS
!      WRITE(6,44771)(FCB(K,IFILEX),K=1,15)
!      WRITE(6,44773)(MEM(IND1+K),K=1,8)
 44773 FORMAT(' DBMMGR,BUFFER,IND1=',8(1X,z8))
 DO  i = 1, nbuff3
   mem( ind2+i ) = mem( ind1+i )
 END DO
 GO TO 7000
 7000  CONTINUE
! SET INDBAS TO POINT TO CURRENT BUFFER
 indbas = fcb( 2, ifilex )
!      IF ( NAME .NE. 307 ) GO TO 7777
!      IF ( IFILEX .NE. 48 ) GO TO 7777
!      PRINT *,' DBMMGR RETURNING,IFILEX,INDBAS=',IFILEX,INDBAS
!      PRINT *,' DBMMGR RETURNING,INDCLR,INDCBP=',INDCLR,INDCBP
!      write(6,40648)(mem(kb),kb=indbas-4,indbas+8)
 40648 FORMAT(' returned buffer=',/,10(4(1X,z8),/))
!      WRITE(6,44771)(FCB(K,IFILEX),K=1,15)
!      CALL DBMFDP
 44771 FORMAT(' returned FCB=',/,2(5I8,/),2I8,4X,2A4,4X,i8)
 7777  CONTINUE
 RETURN
END SUBROUTINE dbmmgr
