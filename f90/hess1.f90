SUBROUTINE hess1 (kdd,mdd,lamd,phid,oeigs,nfound,nvecd,bdd,scr1,  &
        scr2,scr3,scr4,scr5,scr6,scr7,eed,method)
     
!     SUBROUTINE HESS1 TRANSFORMS THE  PROBLEM
!         PSQ M  + P B  + K   INTO   PSQ I  + MINV K
 
!     THREE CASES ARE AVAILABLE
!         1    BDB = 0   MDD  NOT IDENTITY
!              AMAT=     MINVERSE K  (MINUS ADDED IN CORE)
!                        OUTPUT   P = CSQRT COMPUTED  PS
!                        OUTPUT VEC = COMPUTED VECTOR
 
!         2    BDD = 0   MDD  IDENTITY
!              AMAT=     KDD
!                        OUTPUT AS IN CASE 1
 
!         3    BDD NOT  ZERO   MDD  NOT  IDENTITY
!              AMAT=  1    1    1
!                     1  0 1-I  1
!                     1----------
!                     1 -1 1 -1 1
!                     1M  K1M  B1
!                     1    1    1
!                     OUTPUT  P   = COMPUTED  P
!                     OUTPUT  VEC = FIRST HALF OF COMPUTED VECTOR
 
!     CORE  LAYOUT (FOR ALLMAT) IS AS FOLLOWS)
 
!     CONTENTS                SIZE              POINTER   TYPE  NAME
!     --------                ----              -------   ----  ----
!     INPUT MATRIX--VECTORS   2*NROW*NROW        IA       COMP  A
!     EIGENVALUES             2*NROW             IL       COMP  LAMBDA
!     H MATRIX                2*NROW*NROW        IH       COMP  H
!     HL MATRIX               2*NROW*NROW        IHL      COMP  HL
!     VECTOR STORAGE          2*NROW             IV       COMP  VEC
!     MULTPLIERS              2*NROW             IM       COMP  MULT
!     INTH                    NROW               INTH     INT   INTH
!     INT                     NROW               INT      LOG   INT
 
!     BUFFER                  SYSBUF             IBUF1    INT   BUFFER
 
 
!     VARIABLE  DEFINITION
 
!     ID   0  MEANS  IDENTY MASS MATRIX
!     IBDD 0  MEANS  NULL B MATRIX
!     AMAT    FINAL  A MATRIX GINO NAME
!     NROW    ORDER  OF PROBLEM
 
 
 
 
 INTEGER, INTENT(IN)                      :: kdd
 INTEGER, INTENT(IN)                      :: mdd
 INTEGER, INTENT(IN)                      :: lamd
 INTEGER, INTENT(IN)                      :: phid
 INTEGER, INTENT(IN)                      :: oeigs
 INTEGER, INTENT(OUT)                     :: nfound
 INTEGER, INTENT(OUT)                     :: nvecd
 INTEGER, INTENT(IN)                      :: bdd
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN OUT)                  :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER, INTENT(IN OUT)                  :: scr5
 INTEGER, INTENT(IN OUT)                  :: scr6
 INTEGER, INTENT(IN)                      :: scr7
 INTEGER, INTENT(IN)                      :: eed
 INTEGER, INTENT(IN)                      :: method
 INTEGER :: iz(8),mcb(7),sysbuf,NAME(2),FILE,ihead(10), amat,eigc(2),poin
 DOUBLE PRECISION :: d1,d2,d3,d4,d5,dz(1),temp(2)
 COMPLEX :: cz(1),tz
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /unpakx/ itc,ii,jj,incr
 COMMON /zzzzzz/ z(1)
 COMMON /cdcmpx/ dum32(32),ib
 COMMON /system/ ksystm(65)
 COMMON /output/ head(1)
 EQUIVALENCE     (ksystm( 1),sysbuf), (ksystm(2),mout ),  &
     (ksystm(55),iprec ), (z(1),dz(1),iz(1),cz(1))
 DATA    NAME  / 4HHESS,4H1        /
 DATA    ihead / 0,1009,4,7*0      /
 DATA    eigc  , poin/ 207,2,4HPOIN/
 DATA    iz0   / 0   /
 
!     DETERMINE  IF MASS MATRIX IS IDENTITY
 
 mcb(1) = mdd
 CALL rdtrl (mcb)
 id = 0
 IF (mcb(4) == 8) id = 1
 nrow  = mcb(2)
 amat  = kdd
 IF (id /= 0) GO TO 10
 
!     DECOMPOSE  MASS MATRIX
 
 ib = 0
 CALL cfactr (mdd,scr1,scr2,scr3,scr4,scr5,iopt)
 
!     SOLVE FOR AMATRIX
 
 CALL cfbsor (scr1,scr2,kdd,scr3,iopt)
 
!     DETERMINE IF  B MATRIX IS NULL
 
 amat = scr3
 10 ibdd = 0
 mcb(1) = bdd
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0 .OR. mcb(6) == 0) GO TO 30
 
!     FORM  M-1  B
 
 ibdd  = 1
 imat1 = bdd
 imat2 = kdd
 IF (id /= 0) GO TO 20
 
!     - AS OF APRIL 1985 -
!     THE UPPER AND LOWER TRIANGULAR MATRICES IN SCR1 AND SCR2 WERE
!     MYSTERIOUSLY DESTROYED HERE. MUST CALL CFACTR TO RE-GENERATE THEM
 
!     - AS OF JUNE 1991 -
!     TRY WITHOUT 2ND CALL TO CFACTR, AND MAKE SURE SCR1 AND SCR2 ARE
!     STILL GINO UNITS 301 AND 302
 
 ib = 0
 CALL cfactr (mdd,scr1,scr2,scr3,scr4,scr5,iopt)
 
 CALL cfbsor (scr1,scr2,bdd,scr4,iopt)
 imat1 = scr4
 imat2 = scr3
 20 CALL hess2 (nrow,scr5,scr6)
 
!     IDENTITY ON SCR5  MERGE VECTOR ON SCR6
 
 CALL merged (0,scr5,imat2,imat1,scr7,scr6,scr6,0,0)
 amat = scr7
 nrow = 2*nrow
 
!     ALLOCATE  CORE FOR  ALLMAT
 
 30 ia  = 1
 il  = ia + 2*nrow*nrow
 ih  = il + 2*nrow
 ihl = ih + 2*nrow*nrow
 iv  = ihl+ 2*nrow*nrow
 im  = iv + 2*nrow
 inth= im + 2*nrow
 INT = inth + nrow
 nz  =  korsz(iz)
 ibuf1 = nz - sysbuf + 1
 IF (ih+sysbuf > nz) CALL mesage (-8,0,NAME)
 
!     PROCESS EIGC CARD
 
 FILE = eed
 CALL preloc (*900,iz(ibuf1-1),eed)
 CALL locate (*900,iz(ibuf1-1),eigc,iflag)
 50 CALL fread (eed,iz,10,0)
 IF (method == iz(1) .OR. method == -1) GO TO 70
 
!     SKIP REMAINDER OF EIGC CARD
 
 60 CALL fread (eed,iz,7,0)
 IF (iz(6) /= -1) GO TO 60
 GO TO 50
 
!     EIGC  CARD  FOUND
 
 70 inorm = 0
 IF (iz(4)  /= poin) inorm = 1
 isil  = iz(6)
 epsi  = 1.0E-6
 IF (z(iz0+8) /= 0.0) epsi = z(iz0+8)
 
!     PROCESS  REGION  DEFINITION
 
 CALL fread (eed,iz,7,0)
 alph1 = z(1)
 alph2 = z(iz0+3)
 w1    = z(iz0+2)
 w2    = z(iz0+4)
 nvecd = iz(7)
 IF (nvecd > 0) GO TO 95
 
!     ---- SET DEFAULT TO ONE SOLUTION VECTOR ----
 
 nvecd = 1
 WRITE  (mout,90) uwm
 90 FORMAT (a25,' 2357, ONE VECTOR (DEFAULT) WILL BE COMPUTED IN THE',  &
     ' COMPLEX REGION.')
 95 CALL CLOSE (eed,1)
 nvecd = MAX0(nvecd,1)
 
!     BRING IN  TERMS OF MATRIX
 
 CALL gopen (amat,iz(ibuf1),0)
 itc  =-3
 ii   = 1
 jj   = nrow
 incr = 1
 DO  i = ia,il
   z(i) = 0.0
 END DO
 j = ia
 DO  i = 1,nrow
   CALL unpack (*110,amat,z(j))
   110 j = j + 2*nrow
 END DO
 CALL CLOSE (amat,1)
 
!     DO IT
 
 ncount = nvecd
 CALL allmat (z(ia),z(il),z(ih),z(ihl),z(iv),z(im),z(inth),z(INT),  &
     nrow,ncount,iopt1)
 nfound = ncount/iprec
 FILE   = lamd
 CALL OPEN (*900,lamd,iz(ibuf1),1)
 DO  i = 1,nrow
   j  = ia + nrow*nrow + i - 1
   IF (ibdd /= 0) GO TO 210
   
!     PUT OUT COMPLEX SQUARE ROOT
   
   tz = CSQRT(cz(j))
   IF (AIMAG(tz) < 0.0) tz = -tz
   temp(1) = REAL(tz)
   temp(2) = AIMAG(tz)
   GO TO 220
   
!     NON-ZERO  B
   
   210 CONTINUE
   temp(1) = REAL(cz(j))
   temp(2) = AIMAG(cz(j))
   220 CALL WRITE (lamd,temp,4,1)
 END DO
 CALL CLOSE (lamd,1)
 
!     PUT OUT  EIGENVECTORS
 
 FILE = phid
 CALL OPEN (*900,phid,iz(ibuf1),1)
 j    = nrow*nrow + nrow
 k    = ia - 1
 nout = nrow*2
 IF (ibdd /= 0) nout = nout/2
 DO  m = 1,nvecd
   d1 = 0.0
   DO  i = 1,nout,2
     ii = j + i
     jj = k + i
     dz(ii  ) = z(jj  )
     dz(ii+1) = z(jj+1)
     d2 = dz(ii)*dz(ii) + dz(ii+1)*dz(ii+1)
     IF (d2 < d1) CYCLE
     d3 = dz(ii  )
     d4 = dz(ii+1)
     d1 = d2
   END DO
   IF (inorm == 0) GO TO 350
   320 DO  i = 1,nout,2
     jj = j + i
     d5 = (dz(jj)*d3 + dz(jj+1)*d4)/d1
     dz(jj+1) = (d3*dz(jj+1) - d4*dz(jj))/d1
     dz(jj  ) = d5
   END DO
   GO TO 360
   350 jj = 2*isil + j
   d2 = dz(jj)*dz(jj) + dz(jj-1)*dz(jj-1)
   IF (d2 == 0.0D0 .OR. d1/d2 > 1.0D6) GO TO 320
   d3 = dz(jj-1)
   d4 = dz(jj  )
   d1 = d2
   GO TO 320
   360 CONTINUE
   CALL WRITE (phid,dz(j+1),nout*2,1)
   k  = k + nrow*2
 END DO
 CALL CLOSE (phid,1)
 
!     PUT OUT OEIGS
 
 CALL gopen (oeigs,iz(ibuf1),1)
 CALL WRITE (oeigs,ihead,10,0)
 iz(1) = nfound
 iz(2) = nvecd
 iz(3) = 0
 iz(4) = 0
 iz(5) = 0
 iz(6) = 0
 iz(7) = 0
 iz(8) = 1
 CALL WRITE (oeigs,iz,40,0)
 CALL WRITE (oeigs,head,96,1)
 CALL CLOSE (oeigs,1)
 mcb(1) = oeigs
 mcb(2) = nfound
 mcb(3) = nvecd
 CALL wrttrl (mcb)
 RETURN
 
!     ERROR MESSAGES
 
 900 ip1 =-1
 CALL mesage (ip1,FILE,NAME)
 RETURN
 
END SUBROUTINE hess1
