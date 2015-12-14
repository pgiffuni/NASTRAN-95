SUBROUTINE dbmio ( opcode )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 INCLUDE 'XNSTRN.COM'
 
!  OPCODE
!         = 1   OPEN, IOCODE = 0 OPEN FOR READ WITH REWIND
!                            = 1 OPEN FOR WRITE WITH REWIND
!                            = 2 OPEN FOR READ WITHOUT REWIND
!                            = 3 OPEN FOR WRITE WITHOUT REWIND
!         = 2   CLOSE, IOCODE = 1 CLOSE WITH REWIND
!                                OTHERWISE, NO REWIND
!         = 3   REWIND
!         = 4   WRITE ONE BLOCK
!         = 5   READ ONE BLOCK
!         = 6   POSITION
!         = 7   DELETE FILE
!         = 8   WRTBLK CODE
!         = 9   RDBLK CODE
!-----------------------------------------------------------------------
 
 
 INTEGER, INTENT(IN OUT)                  :: opcode
 
!      PRINT *,' DBMIO CALLED WITH OPCODE,IFILEX=',OPCODE,IFILEX
!      PRINT *,' DBMIO,NBLOCK,IOCODE=',NBLOCK,IOCODE
!      WRITE(6,40646)(FCB(K,IFILEX),K=1,15)
 40646 FORMAT(' DBMIO-ENTRY,FCB=',/ i3,i7,4I5,i7,i2,4I7,1X,2A4,i4)
 SELECT CASE ( opcode )
   CASE (    1)
     GO TO 100
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
!-OPEN ------------------------------
!     OPEN FILE ACCORDING TO IOCODE
!       =0, OPEN AND READ FIRST BLOCK
!       =1, OPEN AND RETURN ( OPEN FOR WRITE )
!       =2, OPEN AND READ THE CURRENT BLOCK WHEN FILE WAS CLOSED
!       =3, OPEN AND READ THE CURRENT BLOCK WHEN FILE WAS CLOSED
 100   CONTINUE
 CALL dsgnop
 fcb( 15, ifilex ) = 700+iocode
 fcb(  1, ifilex ) = iocode
 IF ( iocode == 0 ) GO TO 110
 IF ( iocode == 1 ) GO TO 111
 IF ( iocode == 2 ) GO TO 112
 IF ( iocode == 3 ) GO TO 113
 110   CONTINUE
 fcb( 4, ifilex ) = nblock
 GO TO 600
 111   CONTINUE
 fcb( 4, ifilex ) = nblock
 fcb( 5, ifilex ) = nblock
 fcb( 6, ifilex ) = nblock
 GO TO 7000
 112   CONTINUE
 113   CONTINUE
 nblock = fcb( 4, ifilex )
 IF ( fcb( 5, ifilex ) /= 0 ) GO TO 600
 nblock = 1
 GO TO 111
!-CLOSE -----------------------------
 200   CONTINUE
 CALL dsgncl
 IF ( iocode == 0 ) fcb( 4, ifilex ) = 1
 IF ( fcb( 15, ifilex ) /= 701 .AND. fcb( 15, ifilex ) /= 703 ) GO TO 210
 fcb(  4, ifilex ) = fcb(  4, ifilex ) - 1
 fcb(  6, ifilex ) = fcb(  6, ifilex ) - 1
 210   CONTINUE
 fcb( 15, ifilex ) = 0
 GO TO 7000
!-REWIND ----------------------------
 300   CONTINUE
 fcb(  4, ifilex ) = 1
 nblock = 1
 IF ( fcb( 15, ifilex ) == 701 .OR. fcb( 15, ifilex ) == 703 ) GO TO 7000
 GO TO 600
!-WRITE -----------------------------
 400   CONTINUE
 CALL dsgnwr
 fcb( 4, ifilex ) = fcb( 4, ifilex ) + 1
 IF ( fcb( 4, ifilex ) > fcb( 6, ifilex ) )  &
     fcb( 6, ifilex ) = fcb( 4, ifilex )
 GO TO 7000
!-READ
 500   CONTINUE
 fcb( 4, ifilex ) = fcb( 4, ifilex ) + 1
 nblock = fcb( 4, ifilex )
 CALL dsgnrd
 GO TO 7000
!-POSITION AND READ BLOCK "NBLOCK"
 600   CONTINUE
 CALL dsgnrd
 GO TO 7000
!-DELETE FILE
 700   CONTINUE
 OPEN  (ifilex, FILE=mdsnam(ifilex), STATUS='UNKNOWN')
 CLOSE (ifilex, STATUS='DELETE')
 fcb( 5, ifilex ) = 0
 fcb( 6, ifilex ) = 0
 GO TO 7000
!-SPECIAL RDBLK CALL
 900   CONTINUE
 PRINT *,' ERROR, DBMIO CALLED FOR RDBLK CALL'
 STOP
!-SPECIAL WRTBLK CALL
 800   CONTINUE
 PRINT *,' ERROR, DBMIO CALLED FOR WRTBLK CALL'
 STOP
 7000  CONTINUE
!      WRITE(6,40647)(FCB(K,IFILEX),K=1,15)
 40647 FORMAT(' DBMIO-EXIT,FCB=',/ i3,i7,4I5,i7,i2,4I7,1X,2A4,i4)
 RETURN
END SUBROUTINE dbmio
