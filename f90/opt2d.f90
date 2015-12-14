SUBROUTINE opt2d (ipr,pr)
!-----
!   COPY OPTP1 TO OPTP2 DATA FILE.
!  CHANGE RECORD 3      WORD 1 = IABS (PID).
!                       WORD 4 = PLST
!                       WORD 5 = ALPH
!-----
 
 INTEGER, INTENT(OUT)                     :: ipr(1)
 REAL, INTENT(OUT)                        :: pr(1)
 
 INTEGER :: zcor     ,eor      , optp1    ,optp2    ,iz(1)
 
 COMMON /BLANK/ skp1(9),nwdsp,optp1,skp3(2),optp2,skp4(2),nprw
 COMMON /names / nrd,nrrew,nwrt,nwrew,next
 COMMON /optpw2/ zcor,z(1)
 EQUIVALENCE (iz(1),z(1))
 
!  . RECORD ZERO - COPY NAME AND 6 PARAMETERS...
 
 CALL fread (optp1,z(1),8,next)
 CALL fname(optp2,z(1))
 CALL WRITE (optp2,z(1),8,next)
 
!  . RECORD ONE (POINTERS) AND TWO (ELEMENT DATA)...
 
 DO  i = 1,2
   n = zcor
   10 eor = next
   CALL READ(*20,*20,optp1,z,zcor,0,n)
   eor = 0
   20 CALL WRITE (optp2,z(1),n,eor)
   IF (eor == 0) GO TO 10
 END DO
 
!  . RECORD THREE - PROPERTY DATA...
 
 eor = 0
 DO  i = 1,nprw,nwdsp
   ipr(i) = IABS(ipr(i) )
   pr(i+4) = -1.0
   CALL WRITE (optp2,ipr(i),nwdsp,eor)
 END DO
 CALL WRITE (optp2,0,0,next)
 
!  . RECORD FOUR - PLIMIT DATA...
 
 CALL fread (optp1,0,0,next)
 n = zcor
 50 eor = next
 CALL READ(*60,*60,optp1,z,zcor,0,n)
 eor = 0
 60 CALL WRITE (optp2,z(1),n,eor)
 IF (eor == 0) GO TO 50
 
 CALL eof (optp2)
 iz(1) = optp1
 CALL rdtrl(iz(1))
 iz(1) = optp2
 CALL wrttrl (iz(1))
 RETURN
END SUBROUTINE opt2d
