SUBROUTINE rcovsl (NAME,item,in,amat,scr2,scr3,out,z,iz,lcore,  &
        first,rfno)
     
!     RCOVSL CALCULATES THE STATIC LOAD VECTORS FOR THE SUBSTRUCTURING
!     PHASE 2 AND PHASE 3 OPERATIONS FROM THE SUBSTRUCTURE SOLN ITEM
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN)                      :: item
 INTEGER, INTENT(IN)                      :: in
 INTEGER, INTENT(IN)                      :: amat
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(OUT)                     :: out
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(IN)                      :: iz(1)
 INTEGER, INTENT(IN)                      :: lcore
 LOGICAL, INTENT(IN OUT)                  :: first
 INTEGER, INTENT(IN)                      :: rfno
 
 INTEGER :: pmx,fmx,cmx,slmx,t,  &
     signpf,signc,prec,scr,rd,rdrew,wrt,wrtrew,rew,  &
     otypp,sysbuf,soln,srd,subr(2),buf1,fss(2),ibuf(3),  rc, TYPE
 
 COMMON /mpyadx/ pmx(7),fmx(7),cmx(7),slmx(7),mcore,t,signpf,signc, prec,scr
 COMMON /names / rd,rdrew,wrt,wrtrew,rew,norew
 COMMON /packx / itypp,otypp,irowp,nrowp,incp
 COMMON /system/ sysbuf,nout
 DATA    soln  / 4HSOLN /, srd   / 1 /
 DATA    subr  / 4HRCOV ,  4HSL  /
 
!     INITIALIZE
 
 buf1  = lcore - sysbuf + 1
 itypp = 1
 irowp = 1
 incp  = 1
 mcore = lcore
 t     = 0
 signpf= 1
 prec  = 0
 
!     READ LOAD MATRIX FROM SOF ONTO GINO FILE
 
 pmx(1) = in
 CALL rdtrl (pmx)
 IF (pmx(1) > 0) GO TO 5
 itm = item
 CALL mtrxi (scr2,NAME,item,z(buf1),rc)
 IF (rc == 3) GO TO 600
 IF (rc /= 1) GO TO 1000
 pmx(1) = scr2
 CALL rdtrl (pmx)
 5 nrowp = pmx(2)
 TYPE  = pmx(5)
 IF (rfno == 8 .AND. TYPE <= 2) TYPE = TYPE + 2
 otypp = TYPE
 IF (first) GO TO 500
 
!     PROCESS INITIAL SOLN DATA
 
 itm = soln
 CALL sfetch (NAME,soln,srd,rc)
 IF (rc /= 1) GO TO 1000
 CALL suread (fss,2,n,rc)
 IF (rc /= 1) GO TO 1100
 CALL suread (ibuf,3,n,rc)
 IF (rc /= 1) GO TO 1100
 IF (rfno == 3) GO TO 600
 nb  = ibuf(2)
 nst = ibuf(3)
 
!     INTILIZE SCR1 FILE
 
 CALL makmcb (fmx,amat,nrowp,2,TYPE)
 CALL gopen (amat,z(buf1),wrtrew)
 
!     PACK FACTOR MATRIX FOR R. F. 1,2
 
 IF (rfno == 8 .OR. rfno == 9) GO TO 100
 DO  i = 1,nst
   DO  j = 1,nrowp
     z(j) = 0.0
   END DO
   n = 1
   CALL sjump (n)
   IF (n < 0) GO TO 1200
   CALL suread (nl,1,n,rc)
   IF (rc /= 1) GO TO 1100
   IF (nl < 0) CYCLE
   IF (nl == 0) GO TO 30
   IF (nrowp+2*nl >= buf1) CALL mesage (-8,0,subr)
   CALL suread (z(nrowp+1),2*nl,n,rc)
   IF (rc /= 1) GO TO 1100
   nrow = nrowp - 1
   DO  j = 1,nl
     nrow = nrow + 2
     nr   = iz(nrow)
     z(nr)= z(nrow+1)
   END DO
   30 CALL pack (z(1),amat,fmx)
 END DO
 CALL CLOSE (amat,rew)
 CALL wrttrl(fmx)
 GO TO 500
 
!     PACK FACTOR MATRIX FOR R. F. 8,9
 
 100 CALL suread (iz(1),3*nb,n,rc)
 IF (rc /= 1) GO TO 1100
 CALL suread (nl,1,n,rc)
 IF (rc /= 1) GO TO 1100
 IF (nl <= 0) GO TO 600
 IF (nl >= buf1) CALL mesage (-8,0,subr)
 CALL suread (iz(1),nl,n,rc)
 IF (rc /= 1) GO TO 1100
 n  = 1
 CALL sjump (n)
 IF (n < 0) GO TO 1200
 ip = 1
 IF (rfno == 8) ip = 2
 IF (rfno == 8) itypp = 3
 ifact = nl + 1
 nfact = nl + nl*ip
 icol  = nfact + 1
 ncol  = nfact + ip*nrowp
 IF (ncol >= buf1) CALL mesage (-8,0,subr)
 
 DO  i = 1,nst
   DO  j = icol,ncol
     z(j) = 0.0
   END DO
   n = 1
   CALL sjump (n)
   IF (n < 0) GO TO 1200
   CALL suread (z(ifact),nl*ip,n,rc)
   IF (rc /= 1) GO TO 1100
   nrow = ifact - ip
   nrs  = icol  - ip
   DO  j = 1,nl
     nrow = nrow + ip
     nr   = nrs  + iz(j)*ip
     z(nr)= z(nrow)
     IF (ip == 2) z(nr+1) = z(nrow+1)
   END DO
   CALL pack (z(icol),amat,fmx)
 END DO
 CALL CLOSE (amat,rew)
 CALL wrttrl (fmx)
 
!     OUT = LOADS*FACTORS
 
 500 fmx(1) = amat
 CALL rdtrl (fmx)
 cmx(1) = 0
 CALL makmcb (slmx,out,pmx(3),2,TYPE)
 scr = scr3
 CALL mpyad (z,z,z)
 CALL wrttrl (slmx)
 GO TO 700
 
!     NO SCALAR LOADS
 
 600 out = 0
 CALL CLOSE (amat,rew)
 700 RETURN
 
!     ERRORS
 
 1000 CALL smsg (rc-2,itm,NAME)
 GO TO 600
 1100 CALL smsg (rc+4,itm,NAME)
 GO TO 600
 1200 CALL smsg (7,itm,NAME)
 GO TO 600
END SUBROUTINE rcovsl
