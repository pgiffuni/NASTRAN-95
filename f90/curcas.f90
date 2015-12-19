SUBROUTINE curcas(*,nskip,trl,mcb,zz,ibuf)
!   THIS SUBROUTINE COPIES MATRIX FILE TRL(1) TO FILE MCB(1)
! SKIPPING NSKIP-1 MATRIX COLUMNS.  PRIMARY USE IS TO CREATE A MATRIX
! THAT INCLUDES ONLY SUBCASES IN THE CURRENT DMAP LOOP.
!   ALL FILES ARE OPENED, CLOSED AND TRIALERS WRITTEN.
!   IF NSKIP WOULD RESULT IN NO-COPY, MCB(1) IS SET TO TRL(1).
!     TRL - INPUT TRAILER FOR FILE BEING CONVERTED.
!     MCB - OUTPUT TRAILER - WORD 1 HAS GINO FILE NAME.
!     ZZ  - OPEN CORE.
!     IBUF- LOCATION OF TWO GINO BUFFERS.
!     NSKIP - ONE MORE THAN THE SUBCASES TO SKIP.
!       * - NONSTANDARD RETURN IF UNABLE TO PROCESS.
!-----
 
 INTEGER, INTENT(IN)                      :: nskip
 INTEGER, INTENT(IN)                      :: trl(7)
 INTEGER, INTENT(IN OUT)                  :: mcb(7)
 INTEGER, INTENT(IN OUT)                  :: zz(1)
 INTEGER, INTENT(IN)                      :: ibuf
 INTEGER :: parm(4)  , count
 
 COMMON /names / ird,irdrw,iwt,iwtrw, krew,knrw,knerw
 COMMON /system/ isbz
 EQUIVALENCE (icnt,rcnt)
 DATA parm(3),parm(4) / 4HCURC,2HAS  /
 
 parm(2) = trl(1)
 IF (nskip <= 1) GO TO 55
!  . FOR STATICS THE NUMBER OF SUBCASES SKIPPED = NO. COLUMNS SKIPPED.
!  .  OTHER ANALYSIS TYPES NEED TO SUPPLY PROPER VALUE FOR NSKIP...
 i = nskip - 1
 ibf2 = ibuf+isbz
 IF (ibuf <= 0) GO TO 100
 
 CALL rdtrl(trl)
 IF (trl(1) <= 0) GO TO 90
 IF (trl(2) <= i) GO TO 110
 CALL OPEN(*90,trl(1),zz(ibf2),irdrw)
 parm(2) = mcb(1)
 CALL OPEN(*90,mcb(1),zz(ibuf),iwtrw)
 CALL WRITE(mcb(1),mcb(1),2,1)
 parm(2) = trl(1)
 CALL fwdrec(*120,trl(1))
 
 mcb(2) = trl(2) - i
 mcb(3) = trl(3)
 mcb(4) = trl(4)
 mcb(5) = trl(5)
 mcb(6) = trl(6)
 DO  j = 1,i
   CALL fwdrec(*120,trl(1))
 END DO
 CALL cpyfil (trl,mcb,zz,ibuf-1,count)
 rcnt = count
 mcb(7) = icnt
 CALL eof (mcb)
 
 CALL CLOSE (trl(1),krw)
 CALL CLOSE (mcb(1),krw)
 CALL wrttrl (mcb(1))
 GO TO 60
 55 mcb(1) = trl(1)
 60 RETURN
 
!  . ERROR MESSAGES...
 
 90 parm(1) = +1
 GO TO 130
 100 parm(1) = +8
 GO TO 130
 110 parm(1) = +7
 GO TO 130
 120 parm(1) = +2
 
 130 CALL mesage (parm(1),parm(2),parm(3))
 RETURN 1
END SUBROUTINE curcas
