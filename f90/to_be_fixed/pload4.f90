SUBROUTINE pload4 (ibuf5,ido,jopen)
     
!     TO GENERATE PLOAD4 PRESSURE LOAD FOR QUAD4 AND TRIA3 ELEMENTS.
 
!     BOTH ELEMENT TYPES MAY BE PRESENT, OR ONLY ONE OF THE TWO IS
!     PRESENT.
 
!     THIS ROUTINE IS CALLED ONLY BY EXTERN IN SSG1 MODULE, LINK5
 
!     THIS ROUTINE  CALLS PLOD4D OR PLOD4S TO COMPUTE LOAD FOR QUAD4
!     ELEMENTS, AND CALLS T3PL4D OR T3PL4S TO COMPUTE LOAD FOR TRIA3
 
!     IN OVERLAY TREE, THIS ROUTINE SHOULD BE IN PARALLELED WITH FPONT
!     ROUTINE, AND FOLLOWED BY PLOD4D/S AND T3PL4D/S. I.E.
 
!                   ( FPONT
!            EXTERN (        ( PLOD4D  (/ZZSSA1/
!                   ( PLOAD4 ( PLOD4S
!                            ( T3PL4D
!                            ( T3PL4S
 
 
 INTEGER, INTENT(IN)                      :: ibuf5
 INTEGER, INTENT(IN)                      :: ido
 INTEGER, INTENT(OUT)                     :: jopen
 LOGICAL :: allin,debug
 INTEGER :: iz(1),NAME(2),FILE,slt,est,quad4,tria3,t3,q4
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ core(1)
 COMMON /loadx / lcare,slt,idum(5),est
 COMMON /system/ ibuf,nout,jdum(52),iprec
 COMMON /pindex/ iest(45),islt(11)
 COMMON /gpta1 / nelem,last,incr,ielem(1)
 EQUIVALENCE     (core(1),iz(1))
 DATA    quad4 , tria3 , NAME          / 64    , 83    , 4HPLOA,4HD4   /
 DATA    debug / .false. /
 
 
!     T3 AND Q4 KEEP TRACK OF THE PRESENCE OF THE CTRIA3 AND CQUAD4
!     ELEMENTS
 
 t3    = 0
 q4    = 0
 lcore = ibuf5 - ibuf
 ido11 = ido*11
 allin = .false.
 IF (ido11 > lcore) GO TO 400
 IF (debug) WRITE (nout,300)
 300 FORMAT (/,' * PLOAD4 IS CALLED FOR ONE LOAD CASE')
 
!     OPEN CORE IS BIG ENOUGH TO HOLD ALL PLOAD4 DATA.
!     READ THEM ALL INTO CORE
!     (BAD NEWS - OPEN CORE AT THIS TIME IS NOT AVAILABLE)
 
 IF (.NOT.allin) GO TO 400
 
 allin = .true.
 FILE  = slt
 imhere= 350
 CALL READ (*620,*630,slt,core,ido11,0,flag)
 
!     OPEN CORE NOT LARGE ENOUGH TO HOLD ALL PLOAD4 DATA
 
 400 IF (jopen == 1) GO TO 415
 jopen = 1
 FILE  = est
 CALL OPEN  (*610,est,core(ibuf5),0)
 CALL fwdrec (*620,est)
 FILE = est
 410 CALL READ (*430,*560,est,ieltyp,1,0,flag)
 415 IF (ieltyp == quad4) GO TO 440
 IF (ieltyp == tria3) GO TO 445
 420 CALL fwdrec (*430,est)
 GO TO 410
 430 IF (t3+q4 /= 0) GO TO 560
 WRITE  (nout,435) ufm
 435 FORMAT (a23,', PLOAD4 PRESSURE LOAD IS USED WITHOUT THE PRESENCE',  &
     ' OF QUAD4 OR TRIA3 ELEMENT')
 imhere = 435
 GO TO 620
 
 440 IF (q4 >= 1) GO TO 420
 q4 = 1
 IF (debug) WRITE (nout,441) t3
 441 FORMAT (/,'   QUAD4 ELEM FOUND. SETTING Q4 TO 1.  T3 =',i3)
 GO TO 450
 445 IF (t3 == 1) GO TO 420
 t3 = 1
 IF (debug) WRITE (nout,446) q4
 446 FORMAT (/,'   TRIA3 ELEM FOUND. SETTING T3 TO 1.  Q4 =',i3)
 450 j  = incr*(ieltyp-1)
 nwords = ielem(j+12)
 iest(1)= 0
 
 FILE = slt
 ib   = 0
 imhere = 550
 DO  j = 1,ido
   IF (allin) GO TO 460
   jsave = j
   IF (j == 1 .AND. t3+q4 >= 2) GO TO 470
   CALL READ (*620,*630,slt,islt,11,0,flag)
   GO TO 470
   460 DO  i = 1,11
     islt(i) = iz(i+ib)
   END DO
   ib = ib + 11
   470 IF (islt(1)-iest(1) < 0) THEN
     GO TO   550
   ELSE IF (islt(1)-iest(1) == 0) THEN
     GO TO   490
   END IF
   480 CALL READ (*560,*560,est,iest,nwords,0,flag)
   GO TO 470
   
   490 IF (ieltyp == tria3) GO TO 520
   
!     PLOAD4 FOR QUAD4 ELEMENT
   
   IF (debug) WRITE (nout,500) iest(1)
   500 FORMAT (' ==> PROCESS PLOAD4 FOR QUAD ELEM',i8)
   SELECT CASE ( iprec )
     CASE (    1)
       GO TO 505
     CASE (    2)
       GO TO 510
   END SELECT
   505 CALL plod4s
   CYCLE
   510 CALL plod4d
   CYCLE
   
!     PLOAD4 FOR TRIA3 ELEMENT
!     SET ISLT(1) TO NEGATIVE FOR PLOAD4/TRIA3 COMPUTATION
   
   520 IF (debug) WRITE (nout,525) iest(1)
   525 FORMAT (' ==> PROCESS PLOAD4 FOR TRIA3 ELEM',i8)
   islt(1) = -IABS(islt(1))
   SELECT CASE ( iprec )
     CASE (    1)
       GO TO 530
     CASE (    2)
       GO TO 540
   END SELECT
   530 CALL t3pl4s
   CYCLE
   540 CALL t3pl4d
   
 END DO
 
 560 IF (t3+q4 >= 2) GO TO 580
 
!     JUST FINISHED EITHER QUAD4 OR TRIA3 ELEMENT. BACKSPACE EST FILE,
!     AND BACKSPACE SLT FILE IF SLT DATA ARE NOT ALREADY IN CORE.
!     REPEAT PLOAD4 (LOAD TYPE 25) COMPUTAION FOR THE OTHER ELEMENT
!     (TRIA3 OR QUAD4) WHICH WE HAVE NOT YET PROCESSED IN THE FIRST
!     PASS. MUST STEP OVER OTHER LOADS THAT MIGHT BE PRESENT
 
 CALL bckrec (est)
 q4    = q4 + 1
 jsave = 0
 IF (allin) GO TO 410
 
 CALL bckrec (slt)
 imhere = 570
 570 CALL READ (*620,*630,slt,i,1,0,flag)
 IF (i /= 25) GO TO 570
 imhere = 573
 CALL READ (*620,*630,slt,i,1,0,flag)
 IF (i /= ido) GO TO 570
 imhere = 575
 CALL READ (*620,*630,slt,islt,6,0,flag)
 IF (islt(6) /= -1) GO TO 570
 imhere = 577
 CALL READ (*620,*630,slt,islt(7),5,0,flag)
 IF (islt(7) /= 0) GO TO 570
 jsave = 1
 GO TO 410
 
 580 IF (jopen == 1) CALL CLOSE (est,1)
 jopen = 0
 IF (allin .OR. jsave >= ido) GO TO 600
 imhere = 590
 j = (ido-jsave)*11
 CALL READ (*640,*640,slt,0,-j,0,flag)
 600 RETURN
 
 610 j = -1
 GO TO 650
 620 j = -2
 GO TO 650
 630 j = -3
 GO TO 650
 640 j = 1
 650 WRITE  (nout,660) imhere,t3,q4,ido,jsave
 660 FORMAT ('   IMHERE =',i5,'   T3,Q4 =',2I3,'   IDO JSAVE =',2I5)
   CALL mesage (j,FILE,NAME(1))
   GO TO 600
 END SUBROUTINE pload4
