SUBROUTINE flbset
     
!     CONSTRUCTS THE HYDROELASTIC USET VECTOR AND WRITES THE CONECT
!     FILE FOR USE IN CORE ALLOCATION DURING MATRIX ASSEMBLY
 
 EXTERNAL      complf   ,rshift   ,andf     ,orf
 
 LOGICAL :: error
 
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,mpt  &
     ,geom3    ,cstm     ,uset     ,eqexin   ,usetf  &
     ,usets    ,af       ,dkgg     ,fbelm    ,frelm  &
     ,conect   ,afmat    ,afdict   ,kgmat    ,kgdict  &
     ,z        ,group(3) ,two      ,ux       ,uy  &
     ,ufr      ,uz       ,uab      ,ui       ,ua  &
     ,mcb(7)   ,NAME(2)  ,FILE     ,total    ,nam(2)
 
 INTEGER :: andf     ,orf      ,complf   ,rshift
 
!     MACHINE AND HALF WORD
 
 COMMON / machin /       mach     ,ihalf    ,jhalf
 
!     GINO FILES
 
 COMMON / flbfil /       geom2    ,ect      ,bgpdt    ,sil  &
     ,mpt      ,geom3    ,cstm     ,uset ,eqexin   ,usetf    ,usets    ,af  &
     ,dkgg     ,fbelm    ,frelm    ,conect ,afmat    ,afdict   ,kgmat    ,kgdict
 
!     CORE POINTERS
 
 COMMON / flbptr /       error    ,icore    ,lcore    ,ibgpdt  &
     ,nbgpdt   ,isil     ,nsil     ,igrav ,ngrav    ,igrid    ,ngrid    ,ibuf1  &
     ,ibuf2    ,ibuf3    ,ibuf4    ,ibuf5
 
!     OPEN CORE
 
 COMMON / zzzzzz /       z(1)
 
!     POWERS OF TWO
 
 COMMON /two    /        two(32)
 
!     USET PIT POSITIONS
 
 COMMON /bitpos /        bit1(6)  ,ua       ,bit2(16) ,ux  &
     ,uy       ,ufr      ,uz       ,uab ,ui
 
 DATA    NAME   /        4HFLBS,4HET   /
 DATA    nam    /        4HCONE,4HCT   /
 DATA    mcb    /        7*0           /
 
!***********************************************************************
 
!     READ SIL INTO CORE
 
 FILE = sil
 isil = icore
 nz   = igrid - isil - 1
 CALL gopen (sil,z(ibuf1),0)
 CALL READ (*1002,*10,sil,z(isil),nz,0,nsil)
 GO TO 1008
 10 CALL CLOSE (sil,1)
 
!     WRITE OUT CONECT FILE
 
!     FILE 1 - FOR USE IN ASSEMBLING AF MATRIX, CONTAINS SILS WHICH
!              CONNECT FLUID POINTS TO STRUCTURE POINTS ALONG THE
!              BOUNDARY AND SILS WHICH CONNECT FLUID POINTS ALONG THE
!              FREE SURFACE
!     FILE 2 - FOR USE IN ASSEMBLING THE DKGG MATRIX, CONTAINS SILS
!              WHICH CONNECT STRUCTURE POINTS ALONG THE BOUNDARY AND
!              SILS WHICH CONNECT FLUID POINTS ALONG THE FREE SURFACE
 
!     EACH FILE IS COMPOSED OF A 3 WORD RECORD FOR EACH SIL
 
!              WORD      DESCRIPTION
 
!               1        SIL NUMBER
!               2        MAXIMUN GRID POINTS CONNECTED
!               3        MAXIMUM SILS CONNECTED
 
 FILE = conect
 CALL OPEN (*1001,conect,z(ibuf1),1)
 
!     FILE 1
 
 CALL WRITE (conect,nam,2,1)
 DO  i = 1,ngrid
   j = igrid + i - 1
   IF (z(j) <= 0) CYCLE
   nfr = z(j) / 1000000
   nfl = z(j) - nfr*1000000
   group(1) = z(isil+i-1)
   group(2) = nfr+nfl
   group(3) = nfr + 3*nfl
   CALL WRITE (conect,group,3,1)
 END DO
 CALL eof (conect)
 
!     FILE 2
 
 CALL WRITE (conect,nam,2,1)
 DO  i = 1,ngrid
   j = igrid + i - 1
   IF (z(j) >= 0 .AND. z(j) < 1000000) CYCLE
   IF (z(j) > 0) GO TO 30
   ngroup = 3
   nngrid = IABS(z(j))
   nnsil  = nngrid*3
   GO TO 40
   30 ngroup = 1
   nngrid = z(j) / 1000000
   nnsil  = nngrid
   40 jsil   = z(isil+i-1)
   DO  j = 1,ngroup
     group(1) = jsil
     group(2) = nngrid
     group(3) = nnsil
     CALL WRITE (conect,group,3,1)
     jsil = jsil + 1
   END DO
 END DO
 
 CALL CLOSE (conect,1)
 mcb(1) = conect
 mcb(2) = ngrid
 CALL wrttrl (mcb)
 
!     READ USET TABLE INTO CORE
 
 FILE  = uset
 iuset = isil  + nsil  + 1
 nz    = igrid - iuset - 1
 CALL gopen (uset,z(ibuf1),0)
 CALL READ  (*1002,*70,uset,z(iuset),nz,0,nuset)
 GO TO 1008
 70 CALL CLOSE (uset,1)
 
!     CONSTRUCT A LIST OF FREE SURFACE GRID POINTS BY PASSING THROUGH
!     THE GRID POINT CONNECTIVITY TABLE.
 
 icore = iuset + nuset
 ifree = icore
 DO  i = 1,ngrid
   IF (z(igrid+i-1) < 1000000) CYCLE
   z(icore) = i
   icore = icore + 1
   IF (icore >= igrid) GO TO 1008
 END DO
 nfree = icore - ifree
 
!     PASS THROUGH SIL AND PROCESS EACH GRID POINT TO SET THE
!     APPROPRIATE BIT POSITIONS IN THE NEW USET
 
!     *** NOTE.
!     THE UW BIT IS NO LONGER USED.  INSTEAD THE UA BIT WILL REFLECT
!     THE SOLUTION SET  (UAB + UFR)
 
 z(isil+nsil) = nuset + 1
 nstr  = 0
 total = 0
 maska = complf(two(ua))
 juset = iuset
 DO  igrd = 1,nsil
   k = ibgpdt + 4*(igrd-1)
   icstm = z(k)
   IF (icstm == -1) GO TO 82
   IF (z(isil+igrd) == z(isil+igrd-1)+1) GO TO 100
   
!     STURCTURE POINT - SET UX AND UZ. ALSO SET UAB IF UA IS SET
   
   nnsil = 6
   GO TO 84
   82 nnsil = 1
   84 nstr  = nstr + nnsil
   DO  j = 1,nnsil
     z(juset) = orf(z(juset),two(ux))
     z(juset) = orf(z(juset),two(uz))
     IF (andf(z(juset),two(ua)) == 0) GO TO 85
     z(juset) = orf(z(juset),two(uab))
     85 CONTINUE
     total = orf(total,z(juset))
     juset = juset + 1
   END DO
   CYCLE
   
!     FLUID POINT - SET Y BIT.
   
   100 z(juset) = orf(z(juset),two(uy))
   CALL bisloc (*102,igrd,z(ifree),1,nfree,jloc)
   
!     FREE SURFACE FLUID POINT - SET UFR, UA AND UZ BITS
   
   z(juset) = orf(z(juset),two(ufr))
   z(juset) = orf(z(juset),two(ua))
   z(juset) = orf(z(juset),two(uz))
   GO TO 106
   
!     INTERIOR FLUID POINT - SET UI BIT AND TURN OF UA BIT
   
   102 z(juset) = orf(z(juset),two(ui))
   z(juset) = andf(z(juset),maska)
   
   106 CONTINUE
   total = orf(total,z(juset))
   juset = juset + 1
 END DO
 
!     WRITE OUT NEW USETF VECTOR
 
 CALL gopen (usetf,z(ibuf1),1)
 CALL WRITE (usetf,z(iuset),nuset,1)
 CALL CLOSE (usetf,1)
 mcb(1) = usetf
 mcb(2) = 0
 mcb(3) = nuset
 mcb(4) = rshift(total,ihalf)
 mcb(5) = andf(total,jhalf)
 CALL wrttrl (mcb)
 
!     WRITE OUT NEW USETS VECTOR
 
 CALL gopen (usets,z(ibuf1),1)
 luset = iuset + nuset - 1
 DO  i = iuset,luset
   IF (andf(z(i),two(ux)) == 0) CYCLE
   CALL WRITE (usets,z(i),1,0)
 END DO
 CALL CLOSE (usets,1)
 mask   = complf(orf(two(uy),two(ufr)))
 total  = andf(total,mask)
 mcb(1) = usets
 mcb(3) = nstr
 mcb(4) = rshift(total,ihalf)
 mcb(5) = andf(total,jhalf)
 CALL wrttrl (mcb)
 
!     PRINT NEW USET VECTOR IF USER REQUESTS
 
 icore = iuset + nuset
 CALL flbprt (iuset,icore,ibuf1)
 
!     USET PROCESSING COMPLETED
 
 icore = iuset
 RETURN
 
!     ERROR CONDITIONS
 
 1001 n = -1
 GO TO 1100
 1002 n = -2
 GO TO 1100
 1008 n = -8
 1100 CALL mesage (n,FILE,NAME)
 RETURN
END SUBROUTINE flbset
