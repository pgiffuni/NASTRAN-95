SUBROUTINE partn3 (FILE,size,ones,iz,nz,here,buf,core)
!DIR$ INTEGER=64
 
!     CDIR$ IS CRAY COMPILER DIRECTIVE. 64-BIT INTEGER IS USED LOCALLY
!     DO LOOP 10 MAY NOT WORK PROPERLY WITH 48 BIT INTEGER
 
!     PARTN3 CALLED BY PARTN1 AND MERGE1 (VIA PARTN2) BUILDS A BIT
!     STRING AT Z(IZ) THROUGH Z(NZ) AND CONTAINING ONE-BITS ONLY IN
!     THE RESPECTIVE POSITIONS OCCUPIED BY NON-ZERO ELEMENTS IN THE
!     COLUMN VECTOR WHICH IS STORED ON FILE.
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: size
 INTEGER, INTENT(OUT)                     :: ones
 INTEGER, INTENT(IN)                      :: iz
 INTEGER, INTENT(OUT)                     :: nz
 LOGICAL, INTENT(OUT)                     :: here
 INTEGER, INTENT(IN OUT)                  :: buf(4)
 INTEGER, INTENT(IN OUT)                  :: core
 EXTERNAL        lshift,rshift,orf
 LOGICAL :: pass
 INTEGER :: mcb(7),trl(6),bit(64),subr(2),outpt
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /system/ sysbuf,outpt,xxx(37),nbpw
 COMMON /zntpkx/ elem(4),row,eol
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / sym,TYPE,FORM(4),cpcol,rpcol,ireqcl
 EQUIVALENCE     (trl(1),mcb(2))
 DATA    subr  / 4HPART,4HN3  /
 DATA    pass  / .false.      /
 
!     SET UP TABLE OF BITS ON FIRST PASS THROUGH THIS ROUTINE.
 
 IF (pass) GO TO 20
 pass = .true.
 j = nbpw - 1
 k = lshift(1,j)
 DO  i = 1,nbpw
   bit(i) = k
   k = rshift(k,1)
 END DO
 
 20 CALL OPEN (*130,FILE,buf,rdrew)
 here   = .true.
 mcb(1) = FILE
 CALL rdtrl (mcb)
 
!     NUMBER OF WORDS IN COLUMN INCLUDING ZEROS
 
 size = trl(2)
 IF (ireqcl == 0) GO TO 37
 IF (ireqcl > 0 .AND. ireqcl <= trl(1)) GO TO 38
 IF (trl(1) <= 0) GO TO 37
 WRITE  (outpt,30) swm,FILE,trl(1),ireqcl
 30 FORMAT (a27,' 2173, PARTITIONING VECTOR FILE',i5,' CONTAINS',i10,  &
     ' COLUMNS.', /5X,' THE FIRST COLUMN WILL BE USED, NOT THE',  &
     ' REQUESTED COLUMN',i10)
 37 ireqcl = 1
 38 CALL skprec (FILE,ireqcl)
 IF (trl(4) == 1 .OR. trl(4) == 2) GO TO 60
 WRITE  (outpt,50) swm,FILE
 50 FORMAT (a27,' 2174, PARTITIONING VECTOR ON FILE',i5,  &
     ' IS NOT REAL-SINGLE OR REAL-DOUBLE PRECISION.')
 
!     ZERO THE BIT STRING
 
 60 nz = iz + (size-1)/nbpw
 IF (nz > core) CALL mesage (-8,0,subr)
 DO  i = iz,nz
   z(i) = 0
 END DO
 
!     SET UP TO UNPACK THE COLUMN
 
 ones = 0
 eol  = 0
 CALL intpk (*120,FILE,0,1,0)
 GO TO 90
 
!     UNPACK THE ELEMENTS AND TURN ON BITS IN THE BIT STRING.  MAINTAIN
!     COUNT OF BITS IN -ONES-.
 
 80 IF (eol > 0.0) THEN
   GO TO   120
 END IF
 90 CALL zntpki
 IF (row > size) GO TO 100
 k     = row - 1
 zword = k/nbpw + iz
 zbit  = MOD(k,nbpw) + 1
 z(zword) = orf(z(zword),bit(zbit))
 ones  = ones + 1
 GO TO 80
 
!     ELEMENT OF COLUMN LIES OUT OF RANGE INDICATED BY TRAILER
 
 100 WRITE  (outpt,110) sfm,FILE
 110 FORMAT (a25,' 2175, THE ROW POSITION OF AN ELEMENT OF A COLUMN ',  &
     'ON FILE',i5, /5X,'IS GREATER THAN NUMBER OF ROWS ', 'SPECIFIED BY TRAILER.')
 GO TO 160
 
!     BIT STRING IS COMPLETE.
 
 120 CALL CLOSE (FILE,clsrew)
 RETURN
 
!     FILE IS PURGED
 
 130 size = 0
 ones = 0
 here = .false.
 RETURN
 
!     FATAL ERROR
 
 160 CALL CLOSE  (FILE,clsrew)
 CALL mesage (-61,0,subr)
 RETURN
END SUBROUTINE partn3
