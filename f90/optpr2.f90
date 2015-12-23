SUBROUTINE optpr2
     
!     THIS ROUTINE IS THE DRIVER FOR PROPERTY OPTIMIZATION, PHASE 2.
 
!     CALLING SEQUENCE
 
!     OPTPR2  OPTP1,OES1,EST1 / OPTP2,EST2 / V,N,PRINT / V,N,TSTART /
!                                            V,N,COUNT / V,N,CARDNO $
!     WHERE   PRINT  = INPUT/OUTPUT - INTEGER, CALL OFP IF 1, SKIP OFP
!                      IF -1
!             TSTART = INPUT - INTEGER, END TIME AT OPTPR1.
!             COUNT  = INPUT/OUTPUT - INTEGER, ITERATION LOOP COUNTER.
!             CARDNO = INPUT/OUTPUT - INTEGER, PUNCHED CARD COUNT
 
!     LOGICAL         DEBUG
 INTEGER :: PRINT,count,ycor,parm(8),b1,NAME(2),crew,FILE,  &
     sysbuf,outtap,iy(1),NONE(2),optp1,oes1,est1,optp2,  &
     ptrry,est2,b2,ptpty,ptely,ptpry,ptply,tg,tl, tstart,zcor
 REAL :: y(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / PRINT,tstart,count,ncard,skp,ycor,b1,nelop,nwdse,  &
     nwdsp,optp1,oes1,est1,optp2,est2, nelw,nprw,nklw,ntotl,conv
 COMMON /optpw2/ zcor,z(200)
 COMMON /zzzzzz/ core(1)
 COMMON /names / nrd,nrrew,nwrt,nwrew,crew
 COMMON /system/ sysbuf,outtap
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 EQUIVALENCE     (y(1),iy(1),parm(8)), (core(1),MAX,parm(1)),  &
                 (parm(4),iprn ), (parm(7),iprnt)
 DATA    NAME  / 4H opt,4HPR2  /,   NONE / 4H (NO,4HNE)   /
 DATA    ptpty , ptely,ptpry,ptply,ptrry / 5*0 /
!     DATA    DEBUG / .FALSE /
 
 optp1 = 101
 oes1  = 102
 est1  = 103
 optp2 = 201
 est2  = 202
 zcor  = 200
 nwdse = 5
 nwdsp = 6
 
!     LOAD /GPTA1/ ON 1108
 
 CALL delset
 
!     STEP 1.  INITIALIZE AND READ POPT DATA
 
 b1  = korsz(core(1)) - sysbuf + 1
 b2  = b1 - sysbuf
 ycor= b2 -1
 IF (b2 <= 6) GO TO 10
 count = count + 1
 conv  = 0.0
 FILE  = optp1
 CALL OPEN  (*105,FILE,parm(b1),nrrew)
 CALL fread (optp1,parm(1),2,0)
 CALL fread (optp1,parm(1),6,1)
 
!     PARM NOW CONTAINS
 
!       1 = MAX  - MAX NUMBER OF ITERATIONS (I)
!       2 = EPS  - CONVERGENCE TEST (R)
!       3 = GAMA - ITERATION FACTOR (R)
!       4 = IPRN - PRINT CONTROL (I)
!     5,6 = KPUN - PUNCH CONTROL (BCD, YES OR NO)
 
!     NEW PROPERTIES ARE CALCULATED BY,
!     PNEW = (PLST*ALPH) / (ALPH + (1-ALPH)GAMA)
 
!     STEP 2. CHECK TIME TO GO
 
 IF (count > MAX) GO TO 105
 CALL tmtogo (tg)
 IF (tg > 0) GO TO 5
 CALL mesage (45,count,NAME)
 count = 0
 GO TO 110
 5 CALL klock (tl)
 tl = (tl-tstart)/count
 IF (tg <= tl) count = MAX
 iprnt = 0
 
!     STEP 3. READ OPTP1 INTO CORE
 
!     RECORD 1 - POINTERS
 
 ycor = ycor - 7
 IF (ycor < ntypes) GO TO 10
 
!     POINTERS TO OPTIMIZING POINTERS
 
 CALL fread (optp1,y(1),ntypes,0)
 
!     NUMBER OF ELEMENT TYPES THAT MAY BE OPTIMIZED
 
 CALL fread (optp1,nelop,1,0)
 
!     ELEMENT AND PROPERTY POINTERS OF (2,NELOP+1) LENGTH
 
 ycor = ycor - ntypes
 i    = 2*(nelop+1)
 ptpty= ntypes + 1
 IF (ycor < i) GO TO 10
 CALL fread (optp1,y(ptpty),i,1)
 
!     RECORD 2 - ELEMENT DATA
 
 ycor  = ycor  - i
 ptely = ptpty + i
 IF (ycor < nwdse+nwdsp) GO TO 10
 CALL READ (*10,*30,optp1,y(ptely),ycor,1,nelw)
 
!     INSUFFICIENT CORE - PRINT START OF EACH SECTION
 
 
 10 CALL page2 (-3)
 i = ntypes + 1
 WRITE  (outtap,20) ufm,NAME,b1,i,ptpty,ptely,ptpry
 20 FORMAT (a23,' 2289, ',2A4,'INSUFFICIENT CORE (',i10,2H ), /9X,i9,  &
     ' = MATERIAL',i9,' = POINTERS',i9,' = ELEMENTS',i9, ' = PROPERTIES')
 CALL CLOSE (FILE,crew)
 GO TO 100
 
!     RECORD 3 - PROPERTY DATA
 
 30 IF (nelw < nwdse) GO TO 50
 ptpry = ptely + nelw
 ycor  = ycor  - nelw
 IF (ycor < nwdsp) GO TO 10
 CALL READ (*10,*40,optp1,y(ptpry),ycor,1,nprw)
 GO TO 10
 
!     RECORD 4 - PLIMIT DATA
 
 40 IF (nprw < nwdsp) GO TO 50
 ptply = ptpry + nprw
 ycor  = ycor  - nprw
 IF (ycor < 0) GO TO 10
 CALL READ (*10,*70,optp1,y(ptply),ycor,1,nklw)
 GO TO 10
 
!     INSUFFICIENT DATA
 
 50 CALL CLOSE (FILE,crew)
 CALL page2 (-2)
 WRITE  (outtap,60) ufm,NAME
 60 FORMAT (a23,' 2302, SUBROUTINE ',2A4,' HAS NO PROPERTY OR ',  &
     'ELEMENT DATA.')
 GO TO 100
 
!     CLOSE OPTP1 FILE.
!     ALLOCATE AN ARRAY WITH STARTING POINTER PTRRY, OF LENGTH EQUALS TO
!     THE NO. OF PROPERTY CARDS (TO BE USED IN OPT2A, 2B, AND 2C)
!     SET VARIABLE NTOTL TO THE TOTAL LENGTH OF WORDS USED IN OPEN CORE
!     RE-ESTABLISH OPEN CORE UPPER LIMIT, YCOR
 
 70 CALL CLOSE (FILE,crew)
 ptrry = ptply + nklw
 ntotl = ptrry + nprw/nwdsp + 1
 iy(ntotl-1) = -1234567
 IF (ntotl > ycor) GO TO 10
 ycor = b2 - 1
 DO  j = ntotl,ycor
   iy(j) = 0
 END DO
 
!     READ STRESS DATA, SET ALPH
 
 FILE = oes1
 CALL gopen (FILE,parm(b1),nrrew)
 CALL opt2a (iy(ptpty),y(ptely),iy(ptely),y(ptpry),iy(ptpry), y(ptrry))
 IF (iy(ntotl-1) /= -1234567) GO TO 120
 CALL CLOSE (FILE,crew)
 IF (count > MAX) GO TO 105
 
!     SET NEW PROPERTY, CHECK FOR CONVERGENCE
 
 CALL opt2b (iy(ptpry),y(ptpry),y(ptply),y(ptrry))
 
!     CREATE EST2, PUNCH PROPERTIES IF CONVERGED
 
 PRINT = -1
 IF (count >= MAX .OR. count <= 1 .OR. conv == 2.) PRINT = 1
 IF (iprn < 0  .AND. MOD(count,IABS(iprn)) == 0) PRINT = 1
 IF (count > MAX .OR. count < 0) GO TO 90
 IF (count == 1 .OR. count >= MAX .OR. MOD(count,IABS(iprn)) == 0  &
     .OR. conv == 2.) iprnt = 1
 FILE = est1
 CALL OPEN  (*95,FILE,parm(b1),nrrew)
 CALL fread (FILE,NONE(1),2,1)
 FILE = est2
 CALL gopen (FILE,parm(b2),nwrew)
 CALL opt2c (y(ptpty),iy(ptely),iy(ptpry),y(ptpry),y(ptrry))
 CALL CLOSE (FILE,crew)
 CALL CLOSE (est1,crew)
 
!     COPY OPTPR1 TO OPTPR2 - CHANGE RECORD 3
 
 90 IF (count > MAX) GO TO 105
 CALL OPEN (*95,optp1,parm(b1),nrrew)
 FILE = optp2
 CALL OPEN  (*95,FILE,parm(b2),nwrew)
 CALL opt2d (iy(ptpry),y(ptpry))
 CALL CLOSE (FILE,crew)
 CALL CLOSE (optp1,crew)
 GO TO 110
 
!     FILE NOT PRESENT
 
 95 CALL mesage (-1 ,FILE,NAME)
 100 CALL mesage (-61,b2,NAME)
 105 count = -1
 CALL CLOSE (optp1,1)
 110 IF (conv == 2.0) count = MAX
 IF (count <=  0) PRINT = 1
 IF (count ==  0) count =-1
 RETURN
 
 120 WRITE  (outtap,125) ntotl,ptrry
 125 FORMAT (32H0*** rr DIMENSION error/optpr2  ,2I7)
 GO TO 100
END SUBROUTINE optpr2
