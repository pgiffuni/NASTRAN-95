SUBROUTINE amg
     
    !     THIS IS THE MAIN DRIVER FOR AEROELASTIC MATRIX GENERATION
 
    !     NOTES ON NEW METHOD IMPLIMENTATION
    !     1. ACPT FILE WILL BE POSITIONED READY TO READ AN INPUT RECORD
    !        LEAVE FILE READY TO READ NEXT RECORD.
 
    !     2. ALWAYS PACK OUT A COLUMN (REALY A ROW) OF NJ LENGTH
    !        OUTPUT FILE, PACKX, AND TRAILER(MCB) WILL BE SET UP
 
    !     3. YOUR ROW POSITION WILL START AT NROW + 1
 
    !     4. ALWAYS BUMP NROW BY THE NUMBER OF ROWS WHICH EXIST IN
    !        YOUR INPUT RECORD
 
    !     5. COMPUTATIONS FOR AJJK MATRIX WILL HAVE 3 BUFFERS OF CORE USED
    !        COMPUTATIONS FOR OTHER MATRICES WILL HAVE 4 BUFFERS USED
 
    LOGICAL :: debug
    INTEGER :: sysbuf,buf1,buf2,buf3,aero,acpt,ajjl,skj,w1jk,  &
        w2jk,tskj,tw1jk,tw2jk
    DIMENSION       fmach(1),nd(1),NAME(2)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /BLANK / nk,nj
    COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
    COMMON /amgp2 / tw1jk(7),tw2jk(7)
    COMMON /system/ sysbuf,iout
    COMMON /packx / iti,ito,ii,nn,incr
    COMMON /amgbug/ debug
    COMMON /zzzzzz/ iz(1)
    DATA    NAME  / 4HAMG ,4H      /
    DATA    aero  / 101/, acpt /102/, ajjl /201/,  &
        skj   / 202/, w1jk /203/, w2jk /204/
 
    debug =.false.
    CALL sswtch (20,j)
    IF (j == 1) debug =.true.
 
    !     USE IZ TO COMPUTE BUFFERS
 
    icore = korsz(iz)
    ifile = 4*sysbuf + 3*nj
    IF (icore <= ifile) GO TO 460
 
    !     OPEN INPUT STRUCTURAL DATA
 
    icore = icore - sysbuf
    CALL gopen (acpt,iz(icore+1),0)
 
    !     OPEN AND SKIP HEADER ON AERO
 
    ifile = aero
    buf1  = icore - sysbuf
    CALL gopen (aero,iz(buf1+1),0)
 
    !     READ 3 INPUT WORDS INTO COMMON
 
    CALL READ (*450,*450,aero,nd,3,1,n)
 
    !     OPEN OUTPUT FILE FOR AJJL MATRIX, SET UP TRAILER AND WRITE HEADER
 
    buf2  = buf1 - sysbuf
    ifile = ajjl
    CALL OPEN  (*440,ajjl,iz(buf2+1),1)
    CALL fname (ajjl,mcb)
    CALL WRITE (ajjl,mcb,2,0)
    CALL WRITE (ajjl,nj,1,0)
    CALL WRITE (ajjl,nk,1,0)
    buf3  = buf2 - sysbuf
    CALL gopen (skj,iz(buf3+1),1)
    ifile = aero
    CALL READ (*440,*10,aero,iz,buf3,0,n)
    GO TO 460
10  nmk   = n/2
    CALL REWIND (aero)
    CALL fwdrec (*450,aero)
    CALL fwdrec (*450,aero)
    CALL WRITE  (ajjl,nmk,1,0)
    CALL WRITE  (ajjl,iz,n,0)
    ifile = acpt
    iz(1) = 0
    n1    = 2
20  CALL READ (*90,*90,acpt,method,1,0,n)
    iz( 1) = iz(1) + 1
    iz(n1) = method
    SELECT CASE ( method )
        CASE (    1)
            GO TO 30
        CASE (    2)
            GO TO 40
        CASE (    3)
            GO TO 50
        CASE (    4)
            GO TO 50
        CASE (    5)
            GO TO 50
        CASE (    6)
            GO TO 60
        CASE (    7)
            GO TO 70
    END SELECT
 
!     DOUBLET LATTICE METHOD
 
30 CONTINUE
   CALL READ (*440,*450,acpt,mcb,3,1,n)
 
   !     NUMBER OF COLUMNS ADDED EQUAL NUMBER OF BOXES
 
   iz(n1+1) = mcb(3)
   iz(n1+2) = mcb(3)
   GO TO 80
 
   !     DOUBLET LATTICE WITH BODIES
 
40 CALL READ (*440,*450,acpt,mcb,2,1,n)
   iz(n1+1) = mcb(1)
   iz(n1+2) = mcb(2)
   GO TO 80
 
   !     MACH BOX  STRIP THEORY  PISTON THEORY
 
50 CALL READ (*440,*450,acpt,mcb,1,1,n)
   iz(n1+1) = mcb(1)
   iz(n1+2) = mcb(1)
   GO TO 80
 
   !     COMPRESSOR BLADE METHOD
 
60 CALL READ (*440,*450,acpt,mcb,5,1,n)
 
   !     NUMBER OF COLUMNS ADDED IS NJ = NK = (NSTNS*NLINES) FOR THE BLADE
 
   iz(n1+1) = mcb(4)*mcb(5)
   iz(n1+2) = iz(n1+1)
   GO TO 80
 
   !     SWEPT TURBOPROP BLADE METHOD
 
70 CALL READ (*440,*450,acpt,mcb,5,1,n)
 
   !     NUMBER OF COLUMNS ADDED IS NJ = NK = (2*NSTNS*NLINES) FOR THE PROP
 
   iz(n1+1) = 2*mcb(4)*mcb(5)
   iz(n1+2) = iz(n1+1)
80 n1 = n1 + 3
   GO TO 20
90 CALL REWIND (acpt)
   CALL WRITE  (ajjl,iz,n1-1,1)
   mcb(1)  = ajjl
   mcb(2)  = 0
   mcb(3)  = nj
   mcb(4)  = 2
   mcb(5)  = 3
   mcb(6)  = 0
   mcb(7)  = 0
   incr    = 1
   tskj(1) = skj
   tskj(2) = 0
   tskj(3) = nk
   tskj(4) = 2
   tskj(5) = 3
   tskj(6) = 0
   tskj(7) = 0
   ifile   = acpt
 
   !     READ MACH NUMBER AND REDUCED FREQUENCY AND LOOP UNTIL COMPLETED
 
100 CALL READ (*210,*210,aero,fmach,2,0,n)
 
   !     NUMBER OF ROWS ADDED BY EACH RECORD ON ACPT
 
   nrow = 0
   isk  = 1
   nsk  = 0
 
   !     SKIP HEADER
 
   CALL fwdrec (*450,acpt)
 
   !     READ A RECORD AND LOOP BY METHOD UNTIL EOF
   !     NSK IS BUMPED BY DRIVERS = COLUMNS BUILT  ISK = NEXT COLUMN
 
110 CALL READ (*200,*200,acpt,method,1,0,n)
   SELECT CASE ( method )
       CASE (    1)
           GO TO 120
       CASE (    2)
           GO TO 130
       CASE (    3)
           GO TO 140
       CASE (    4)
           GO TO 150
       CASE (    5)
           GO TO 160
       CASE (    6)
           GO TO 170
       CASE (    7)
           GO TO 180
   END SELECT
 
   !     DOUBLET LATTICE METHOD
 
120 CALL dlamg (acpt,ajjl,skj)
   GO TO 190
 
   !     DOUBLET LATTICE WITH BODIES
 
130 CALL dlamby (acpt,ajjl,skj)
   GO TO 190
 
   !     MACH BOX
 
140 CALL mbamg (acpt,ajjl,skj)
   GO TO 190
 
   !     STRIP THEORY
 
150 CALL stpda (acpt,ajjl,skj)
   GO TO 190
 
   !     PISTON THEORY
 
160 CALL pstamg (acpt,ajjl,skj)
   GO TO 190
 
   !     COMPRESSOR BLADE METHOD
 
170 CALL amgb1 (acpt,ajjl,skj)
   GO TO 190
 
   !     SWEPT TURBOPROP BLADE METHOD
 
180 CALL amgt1 (acpt,ajjl,skj)
190 IF (nsk > nk) GO TO 400
   IF (nrow > nj) GO TO 420
   GO TO 110
200 CALL REWIND (acpt)
   GO TO 100
210 CALL CLOSE  (aero,1)
   CALL CLOSE  (ajjl,1)
   CALL CLOSE  (skj,1)
   CALL wrttrl (tskj)
   CALL wrttrl (mcb)
 
   !     COMPUTE W1JK - W2JK
 
 
   !     OPEN OUTPUT FILES
 
   CALL fwdrec (*450,acpt)
   ifile = w1jk
   CALL gopen (w1jk,iz(buf1+1),1)
   ifile = w2jk
   CALL gopen (w2jk,iz(buf2+1),1)
   ifile = acpt
 
   !     SET UP PACKX AND TRAILERS
 
   incr = 1
   iti  = 1
   ito  = 1
 
   !     II AND NN ARE BUMPED BY METHOD DRIVERS
 
   ii   = 1
   DO  i = 2,7
       tw1jk(i) = 0
       tw2jk(i) = 0
   END DO
   tw1jk(1) = w1jk
   tw2jk(1) = w2jk
   tw1jk(3) = nk
   tw1jk(4) = 2
   tw1jk(5) = 1
   tw2jk(3) = nk
   tw2jk(4) = 2
   tw2jk(5) = 1
 
   !     READ A RECORD AND LOOP ON METHOD UNTIL EOR
 
230 CALL READ (*300,*300,acpt,method,1,0,n)
   SELECT CASE ( method )
       CASE (    1)
           GO TO 240
       CASE (    2)
           GO TO 250
       CASE (    3)
           GO TO 260
       CASE (    4)
           GO TO 260
       CASE (    5)
           GO TO 260
       CASE (    6)
           GO TO 270
       CASE (    7)
           GO TO 280
   END SELECT
 
   !     DOUBLET LATTICE METHOD
 
240 CALL dlpt2 (acpt,w1jk,w2jk)
   GO TO 290
 
   !     DOUBLET LATTICE WITH BODIES
 
250 CALL dlbpt2 (acpt,w1jk,w2jk)
   GO TO 290
 
   !     STRIP THEORY     PISTON THEORY
   !     MACH BOX
 
260 CALL stppt2 (acpt,w1jk,w2jk)
   GO TO 290
 
   !     COMPRESSOR BLADE METHOD
 
270 CALL amgb2 (acpt,w1jk,w2jk)
   GO TO 290
 
   !     SWEPT TURBOPROP BLADE METHOD
 
280 CALL amgt2 (acpt,w1jk,w2jk)
290 IF (nn > nk) GO TO 410
   GO TO 230
 
   !     DONE
 
300 CALL CLOSE  (acpt,1)
   CALL CLOSE  (w1jk,1)
   CALL CLOSE  (w2jk,1)
   CALL wrttrl (tw1jk)
   CALL wrttrl (tw2jk)
   RETURN
 
   !     ERROR MESSAGES
 
   !     NROW IN RECORDS DID NOT MATCH NJ PARAMETER
 
400 nrow = nsk
   nj   = nk
   GO TO 420
410 nrow = nn
   nj   = nk
420 WRITE  (iout,430) sfm,nrow,nj
430 FORMAT (a25,' 2264, NUMBER OF ROWS COMPUTED (',i4,') WAS GREATER',  &
       ' THAN SIZE REQUESTED FOR OUTPUT MATRIX (',i4,2H).)
   CALL mesage (-61,n,NAME)
440 nms = -1
   GO TO 470
450 nms = -2
   GO TO 470
460 nms = -8
470 CALL mesage (nms,ifile,NAME)
   RETURN
END SUBROUTINE amg
