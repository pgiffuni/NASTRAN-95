SUBROUTINE mintrp(ni,xi,nd,xd,TYPE,symm1,symk1,dz,infile,outfil,  &
        scr,scr1,g,ncore,nogo,ipres)
     
 
 INTEGER, INTENT(IN)                      :: ni
 REAL, INTENT(OUT)                        :: xi(1)
 INTEGER, INTENT(IN)                      :: nd
 REAL, INTENT(OUT)                        :: xd(2)
 INTEGER, INTENT(IN OUT)                  :: TYPE
 INTEGER, INTENT(OUT)                     :: symm1
 INTEGER, INTENT(OUT)                     :: symk1
 REAL, INTENT(IN OUT)                     :: dz
 INTEGER, INTENT(IN)                      :: infile
 INTEGER, INTENT(IN OUT)                  :: outfil
 INTEGER, INTENT(IN OUT)                  :: scr
 INTEGER, INTENT(IN)                      :: scr1
 REAL, INTENT(IN OUT)                     :: g(1)
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER, INTENT(IN)                      :: ipres
 INTEGER :: sysbuf, isng, scrm
 INTEGER :: a,b,c,d,NAME(2), buff,gpoint
 LOGICAL :: nimag
 LOGICAL :: spec
 COMPLEX :: alpha
 DOUBLE PRECISION :: ar,ai
 
 
 COMMON /system/ sysbuf
 COMMON / packx/ iti,ito,ii,nn,incr
 COMMON /mpyadx/ a(7),b(7),c(7),d(7),nwords,nt,isab,isc,ipre,scrm
 COMMON /saddx / nmat,lcore,ma(7),ita,alpha(2),dum(48),mc(7)
 COMMON /unpakx/ iout,in,nnn,incru
 
 EQUIVALENCE (alpha(1),ar),(alpha(2),ai)
 
 DATA NAME /4HMINT,4HRP  /
 
!-----------------------------------------------------------------------
 
 spec = .false.
 nogo = 0
 
!     DETERMINE TYPE OF CALL FOR G
!     NEGATIVE VALUE FOR KO CALL LSPLIN, POSITIVE CALL SSPLIN
 
 
!     CHECK CORE NEED AT LEAST 1 BUFFER + G
 
 ity = IABS(TYPE)
 kd = 0
 IF(ity > 3) kd=1
 ncol = (1+kd)*nd
 IF(sysbuf+ncol*ni  > ncore) CALL mesage(-8,0,NAME)
 
!     PROTECT AGAINST BAD CALL
 IF(symk1 < 0) symk1 = -1
 IF(symm1 < 0) symm1 = -1
 IF(symk1 > 0) symk1 = 1
 IF(symm1 > 0) symm1 = 1
!     TRANSPOSE FLAG ON
 kt = 1
!     SPECIAL CASE
 IF(nd == 1.AND.ity < 4) GO TO 300
 10 IF(TYPE < 0) GO TO 100
 CALL ssplin(ni,xi,nd,xd,symm1,symk1,kd,kt,dz,g,ncore,isng)
 IF(isng == 2) GO TO 999
 GO TO 200
 100 nii = 2*ni
 DO  i=1,nii,2
   xi(i) = 0.0
 END DO
 nii = 2*nd
 DO  i = 1,nii,2
   xd(i) = 0.0
 END DO
 CALL lsplin(ni,xi,nd,xd,symk1,kd,kt,dz,-1.0,-1.0,1.0,g,ncore,isng)
 IF(isng == 2) GO TO 999
!     PUT OUT G
 200 buff = ncore-sysbuf+1
 nimag = .true.
 IF(ity == 3.OR.ity == 6) nimag = .false.
 IF(nimag) GO TO 210
 iti = scr
 scr = outfil
 outfil = iti
 210 ito = 1
 jj = ncol
 iti = 1
 nn = ni
 b(3) = ni
 b(5) = 1
 gpoint = 1
 215 incr = 1
 j = 1
 ii = 1
 b(1) = scr
 b(2) = 0
 b(4) = 2
 b(6) = 0
 b(7) = 0
 CALL gopen(scr,g(buff),1)
 DO  i = j,jj
   CALL pack(g(gpoint),scr,b)
   gpoint = gpoint + ni
 END DO
 CALL CLOSE(scr,1)
 CALL wrttrl(b)
 IF(spec) GO TO 1000
 
!     MULT INFILE BY G
 
 c(1) = 0
 a(1) = infile
 CALL rdtrl(a)
 d(1) = outfil
 d(3) = a(3)
 d(4) = 2
 d(5) = a(5)
 IF(ity == 2.OR.ity == 5) d(5) = 1
 IF(d(5) == 1.AND.a(5) == 4) d(5) = 2
 nwords = ncore
 nt = 0
 isab = 1
 ipre = ipres
 scrm = scr1
 CALL mpyad(g,g,g)
 CALL wrttrl(d)
 IF(nimag) GO TO 1000
 
!     IMAG PART ONLY WANTED
 
 nmat = 1
 lcore = ncore
 ma(1) = outfil
 CALL rdtrl(ma)
 ita = 3
 alpha(1) = (0.0,-1.0)
 mc(1) = scr
 mc(2) = ma(2)
 mc(3) = ma(3)
 mc(4) = 2
 mc(5) = ma(5)
 mc(6) = 0
 mc(7) = 0
 ai = -1.0D0
 IF(ma(5) == 4) ita = 4
 IF(ita == 4) ar = 0.0D0
 CALL sadd(g,g)
 CALL wrttrl(mc)
 GO TO 1000
 
!     TEST FOR SPECIAL CASE
 
 300 nii = 2*ni
 k = 0
 DO  i = 1,nii,2
   k = k+1
   IF(xi(i) == xd(1).AND.xi(i+1) == xd(2)) GO TO 315
 END DO
 GO TO 10
 
!     PACK OUT COLUMN OF INFILE
 
 315 a(1) = infile
 CALL rdtrl(a)
 buff = ncore-sysbuf +1
 CALL gopen(infile,g(buff),0)
 incru = 1
 in = 1
 nnn = a(3)
 iout = a(5)
 IF(k == 1) GO TO 330
 k = k-1
 CALL skprec(infile,k)
 330 CALL unpack(*998,infile,g)
 CALL CLOSE(infile,1)
 spec = .true.
 scr = outfil
 iti = a(5)
 nn = a(3)
 jj = 1
 gpoint = 1
 IF(ity == 3) gpoint = 2
 ito = 1
 IF(ity == 1) ito = 3
 IF(a(5) == 4) ito = ito+1
 b(3) = a(3)
 b(5) = ito
 GO TO 215
 998 CALL mesage(-7,0,NAME)
 999 nogo = 1
 1000 RETURN
END SUBROUTINE mintrp
