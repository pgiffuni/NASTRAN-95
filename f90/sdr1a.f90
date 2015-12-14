SUBROUTINE sdr1a (INPUT,iout)
     
!     THIS ROUTINE MAKES PS AND IUF COMPATABLE TO COMPUTE QS IN
!     CASE OF TRANSIENT ANALYSIS
 
 
 INTEGER, INTENT(IN)                      :: INPUT
 INTEGER, INTENT(IN)                      :: iout
 INTEGER :: sysbuf,bcd1(2),mcb(7),ps,ia(7),core(1100)
 
 COMMON /zzzzzz/ corex(1)
 COMMON /system/ sysbuf,ksystm(65)
 COMMON /BLANK / loadnn
 COMMON /unpakx/ it1,ii,jj,incr
 COMMON /packx / it2,it3,ii1,jj1,incr1
 EQUIVALENCE     (core(1),corex(1))
 
 DATA    bcd1  / 4HSDR1,4HA    /
 
 nz = korsz(core) - sysbuf
 CALL OPEN (*40,INPUT,core(nz+1),0)
 CALL skprec (INPUT,1)
 nz = nz - sysbuf
 loadnn = MAX0(loadnn,1)
 IF (loadnn == 1) GO TO 50
 ia(1) = iout
 CALL rdtrl (ia)
 IF (ia(2) == 0) GO TO 50
 ia(1) = INPUT
 CALL rdtrl(ia)
 IF(ia(2) == 0) CALL mesage (-7,0,bcd1)
 
!     POSITION TO END
 
 CALL gopen  (iout,core(nz+1),0)
 CALL skpfil (iout,+1)
 CALL skpfil (iout,-1)
 CALL CLOSE  (iout,+2)
 ia(1) = iout
 CALL rdtrl (ia)
 ia(7) = 0
 CALL gopen (iout,core(nz+1),3)
 10 mcb(1) = INPUT
 CALL rdtrl (mcb)
 k   = mcb(2)
 it1 = mcb(5)
 it2 = it1
 it3 = it2
 incr  = 1
 incr1 = 1
 DO  i = 1,k
   ii = 0
   CALL unpack (*20,INPUT,core)
   ii1 = ii
   jj1 = jj
   CALL pack (core,iout,ia)
   CYCLE
   20 ii1 = 1
   jj1 = 1
   core(1) = 0
   core(2) = 0
   core(3) = 0
   core(4) = 0
   CALL pack (core,iout,ia)
 END DO
 CALL CLOSE (INPUT,1)
 CALL CLOSE (iout,1)
 CALL wrttrl (ia)
 40 RETURN
 
!     FIRST TIME
 
 50 CALL gopen (iout,core(nz+1),1)
 ia(1) = INPUT
 CALL rdtrl (ia)
 ia(2) = 0
 ia(6) = 0
 ia(7) = 0
 ia(1) = iout
 GO TO 10
 
!     SDR1D -
 
 ENTRY sdr1d (ps,iuf,iuf1,itran)
!     ===============================
 
 IF (itran == 0) GO TO 60
 itran  = 1
 mcb(1) = ps
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 100
 ncolps = mcb(2)
 mcb(1) = iuf
 CALL rdtrl (mcb)
 IF (ncolps == mcb(2)) RETURN
 
!     THIS IS A TRANSIENT PROBLEM
 itran = 0
 
 60 mcb(1) = iuf
 CALL rdtrl (mcb)
 ncolps = mcb(2)/3
 ibf = korsz(core) - sysbuf
 CALL gopen (iuf,core(ibf),0)
 ibf1 = ibf - sysbuf
 CALL gopen (iuf1,core(ibf1),1)
 it1  = mcb(5)
 it2  = it1
 it3  = it2
 incr = 1
 incr1  = 1
 mcb(1) = iuf1
 mcb(2) = 0
 mcb(6) = 0
 mcb(7) = 0
 DO  i = 1,ncolps
   ii  = 0
   CALL unpack (*70,iuf,core)
   ii1 = ii
   jj1 = jj
   GO TO 80
   70 core(1) = 0
   core(2) = 0
   core(3) = 0
   core(4) = 0
   ii1 = 1
   jj1 = 1
   80 CALL skprec (iuf,2)
   CALL pack (core,iuf1,mcb)
 END DO
 CALL CLOSE (iuf1,1)
 CALL CLOSE (iuf,1)
 CALL wrttrl (mcb)
 100 RETURN
END SUBROUTINE sdr1a
