SUBROUTINE cead
     
    !     COMPLEX  EIGENVALUE EXTRACTION  MODULE
 
    !     5  INPUT  FILES -  KDD,BDD,MDD,EED,CASECC
    !     4  OUTPUT FILES -  PHID,LAMD,OCEIGS,PHIDL
    !     12 SCRATCHES FILES
    !     1  PARAMETER
 
    IMPLICIT INTEGER (a-z)
    REAL :: eps
    DIMENSION       eigc(2),error(2),NAME(2),mcb(7),kz(1)
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm,uwm,uim
    COMMON /cinvpx/ ik(7),im(7),ib(7),ilam(7),iphi(7),  &
        idmpfl,iscr(11),noreg,eps,reg(7,10),phidli
    COMMON /BLANK / nfound
    COMMON /system/ sysbuf,nout
    COMMON /zzzzzz/ iz(1)
    EQUIVALENCE     (kz(1),iz(1))
    DATA    NAME  / 4HCEAD,4H    /
    DATA    hes   / 4HHESS/
    DATA    feer  / 4HFEER/
    DATA    error / 4HEED ,4HCEAD/
    DATA    kdd   , bdd,mdd,eed,casecc / 101   , 102,103,104,105    /
    DATA    phid  , lamd,oceigs,phidl  / 201   , 202, 203,   204    /
    DATA    scr1  , scr2,scr3,scr4,scr5,scr6,scr7,scr8,scr9 /  &
        301   , 302, 303 ,304 ,305 ,306 ,307 ,308 ,309  /
    DATA    scr10 , scr11,scr12 / 310   , 311,  312   /
    DATA    det   , inv,eigc(1),eigc(2) /4HDET ,4HINV ,207,2/
    DATA    iz2   , iz6,iz148   /2,6,148/
 
    !     FIND SELECTED EIGC CARD IN CASECC
 
    ibuf  = korsz(iz) - sysbuf
    CALL OPEN (*1,casecc,iz(ibuf),0)
    CALL skprec (casecc,1)
    CALL fread (casecc,iz,166,1)
    CALL CLOSE (casecc,1)
    j = 148
    method = iz(j)
    scr10  = 310
    GO TO 2
1   method = -1
2   FILE = eed
    CALL preloc (*90,iz(ibuf),eed)
    CALL locate (*130,iz(ibuf),eigc(1),iflag)
10  CALL READ (*110,*140,eed,iz(1),10,0,iflag)
    IF (method == iz(1) .OR. method == -1) GO TO 30
20  CALL fread (eed,iz,7,0)
    j = 6
    IF (iz(j) /= -1) GO TO 20
    GO TO 10
 
    !     FOUND DESIRED  EIGC CARD
 
30  CALL CLOSE (eed,1)
    j = 2
    capp = iz(j)
    IF (capp ==  det) GO TO 50
    IF (capp ==  inv) GO TO 40
    IF (capp ==  hes) GO TO 52
    IF (capp == feer) GO TO 45
    GO TO 130
 
    !     INVERSE POWER--
 
40  ik(1) = kdd
    CALL CLOSE (eed,1)
    CALL rdtrl (ik)
    im(1) = mdd
    CALL rdtrl (im)
    ib(1) = bdd
    CALL rdtrl (ib)
    IF (ib(1) < 0) ib(1) = 0
    IF (ib(6) == 0) ib(1) = 0
    ilam(1)  = scr8
    iphi(1)  = scr9
    idmpfl   = oceigs
    iscr( 1) = scr1
    iscr( 2) = scr2
    iscr( 3) = scr3
    iscr( 4) = scr4
    iscr( 5) = scr5
    iscr( 6) = scr6
    iscr( 7) = scr7
    iscr( 8) = lamd
    iscr( 9) = phid
    iscr(10) = scr10
    iscr(11) = scr11
    phidli   = scr12
    eps      = .0001
    CALL cinvpr (eed,method,nfound)
    nvect = nfound
    GO TO 60
 
!     FEER METHOD
 
45 CONTINUE
   CALL cfeer (eed,method,nfound)
   nvect = nfound
   GO TO 60
 
   !     DETERMINANT
 
50 CALL cdetm (method,eed,mdd,bdd,kdd,scr8,scr9,oceigs,nfound,scr1,  &
       scr2,scr3,scr4,scr5,scr6,scr7,scr10)
   nvect = nfound
   GO TO 60
 
   !     HESSENBURG METHOD
 
52 CONTINUE
   mcb(1) = kdd
   CALL rdtrl (mcb)
   nrow   = mcb(2)
   mcb(1) = bdd
   CALL  rdtrl (mcb)
   IF (mcb(1) > 0) nrow = nrow*2
   nz = korsz(kz)
 
   !     IF INSUFFICIENT CORE EXISTS FOR HESSENBURG METHOD.  DEFAULT TO
   !     INVERSE POWER.
 
   IF (6*nrow*nrow+nrow*8 <= nz) GO TO 55
   WRITE  (nout,53) uim
53 FORMAT (a29,' 2365, INSUFFICIENT CORE EXISTS FOR HESSENBURG ',  &
       'METHOD.  CHANGING TO INVERSE POWER OR FEER.')
   GO TO 40
 
   !     SUFFICIENT CORE.  PROCEED WITH HESSENBURG METHOD
 
55 CONTINUE
   CALL hess1 (kdd,mdd,scr8,scr9,oceigs,nfound,nvect,bdd,scr1,scr2,  &
       scr3,scr4,scr5,scr6,scr7,eed,method)
   nfound = nvect
 
   !     LAMD ON SCR8, PHID ON SCR9
 
   !     SORT EIGENVALUES AND PREPARE OUTPUT FILES
 
60 IF (nfound /= 0) GO TO 70
   nfound = -1
   GO TO 80
70 CALL cead1a (scr8,scr9,phidli,lamd,phid,phidl,nfound,nvect,capp)
80 RETURN
 
   !     ERROR MESAGES
 
90 ip1 = -1
100 CALL mesage (ip1,FILE,NAME)
110 ip1 = -2
   GO TO 100
130 ip1 = -7
   GO TO 100
140 CALL mesage (-31,method,error(1))
   GO TO 140

END SUBROUTINE cead
