SUBROUTINE frd2i (fl,nfreq,ncore,qhhl,scr2,scr1,scr3,scr4,nrow)
     
 
 REAL, INTENT(OUT)                        :: fl(1)
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(IN)                      :: qhhl
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER, INTENT(IN OUT)                  :: scr4
 INTEGER, INTENT(OUT)                     :: nrow
 INTEGER :: trl(7),out
 DIMENSION  mcb(7),NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / bov,q,rm
 COMMON /condas/ pi,twopi
 COMMON /system/ isys,out,dum(52),iprec
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /TYPE  / p(2),iwc(4)
 DATA    NAME  / 4HFRD2,4HI   /
 DATA    nhfrdi/ 4HFRDI/
 
 ibuf1 = ncore - isys
 ibuf2 = ibuf1 - isys
 nrow  = 0
 incr  = 1
 incr1 = 1
 ii    = 1
 inn   = 1
 mcb(1)= qhhl
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 1000
 nrow  = mcb(3)
 ni    =(mcb(2)/mcb(3))*2
 nnn   = nrow
 nn    = nrow*nrow
 iti   = 3
 ito   = iti
 iout  = iti
 nwc   = iwc(iti)
 iscr  = scr1
 nloop = 1
 indx  = 0
 xm    = rm
 IF (rm >= 0.0) GO TO 5
 iscr  = scr2
 nloop = nfreq
 indx  = 1
 5 CALL makmcb (trl,iscr,nn,mcb(4),ito)
 
!     MAKE INDEPENDENT FREQ LIST
 
 ipd   = 1
 nl    = 2*nfreq
 n     = nfreq + 1
 icore = ibuf1
 ipi   = ipd + nl
 DO  i = 1,nfreq
   fl(nl) = fl(n-i)*twopi*bov
   fl(nl-1) = 0.0
   nl = nl -2
 END DO
 
!     MAKE INDEPENDENT FREQ LIST
 
 CALL OPEN  (*1000,qhhl,fl(ibuf2),0)
 CALL gopen (iscr,fl(ibuf1),1)
 CALL READ  (*999,*999,qhhl,fl(ipi),-3,0,flag)
 CALL READ  (*999,*999,qhhl,n,1,0,flag)
 n = n + n
 IF (rm >= 0.0 .OR. n == ni) GO TO 15
 WRITE  (out, 2000) ufm,n,ni
 2000 FORMAT (a23,', THE NUMBER OF (M,K) PAIRS SPECIFIED ON MKAEROX ',  &
     'CARDS (', i5, ') IS NOT EQUAL ', /5X,  &
     'TO THE NUMBER OF FREQUENCIES SPECIFIED (', i5, '),')
 CALL mesage (-37,0,NAME)
 15  ni = MIN0(ni,n)
 CALL READ (*999,*999,qhhl,fl(ipi),ni,1,flag)
 IF (rm < 0.0) CALL CLOSE (qhhl, 1)
 
 DO  kkk = 1, nloop
   IF (rm >= 0.0) GO TO 20
   xm = fl(2*kkk)
   CALL gopen (qhhl,fl(ibuf2),0)
   20 CONTINUE
   
!     FOR RM.GE.0.0, FIND M CLOSEST TO XM
!     FOR RM.LT.0.0, FIND K CLOSEST TO XM
   
   icp = ipi + ni
   rmi = 1.e20
   rms = 0.0
   DO  i = 1,ni,2
     rmx = ABS(fl(ipi+i+indx-1)-xm)
     rmi = AMIN1(rmi,rmx)
     IF (rmx > rmi) CYCLE
     rms = fl(ipi+i+indx-1)
   END DO
   rmi = rms
   
!     FOR RM.GE.0.0, SELECT ALL K'S ASSOCIATED WITH RMI
!     FOR RM.LT.0.0, SELECT THE K EQUAL TO RMI
   
   k = 0
   DO  i = 1,ni,2
     IF (fl(ipi+i+indx-1) == rmi) GO TO 120
     
!     SKIP MATRIX
     
     CALL skprec (qhhl,nrow)
     CYCLE
     
!     MAKE MATRIX INTO COLUMN
     
     120 fl(ipi+k+1) = fl(ipi+i)
     k  = k + 2
     ji = icp
     n  = nrow*nwc
     DO  j = 1,nrow
       CALL unpack (*131,qhhl,fl(ji))
       GO TO 135
       131 CALL zeroc (fl(ji),n)
       135 ji = ji + n
     END DO
     
!     DIVIDE IMAG PART OF QHHL BY FREQUENCY
     
     jj = icp + 1
     kk = ji  - 1
     DO  j = jj,kk,2
       fl(j) = fl(j)/fl(ipi+i)
     END DO
     IF (rm < 0.0) fl(ipi+i) = -10000.0
     CALL pack (fl(icp),iscr,trl)
     IF (rm < 0.0) EXIT
   END DO
   150 CALL CLOSE (qhhl,1)
   CALL CLOSE (iscr,1)
 END DO
 
 CALL wrttrl (trl)
 CALL bug (nhfrdi,200,k ,1)
 CALL bug (nhfrdi,200,nfreq,1)
 CALL bug (nhfrdi,200,fl(1),icp)
 IF (rm < 0.0) RETURN
 
!     SETUP TO CALL MINTRP
 
 ni   = k/2
 nogo = 0
 nc   = ncore - icp
 CALL dmpfil (-scr1,fl(icp),nc)
 im   = 0
 ik   = 1
 CALL mintrp (ni,fl(ipi),nfreq,fl(ipd),-1,im,ik,0.0,scr1,scr2,scr3,  &
     scr4,fl(icp),nc,nogo,iprec)
 IF (nogo == 1) GO TO 998
 CALL dmpfil (-scr2,fl(icp),nc)
 RETURN
 
 998 WRITE  (out,9980) ufm
 9980 FORMAT (a23,' 2271, INTERPOLATION MATRIX IS SINGULAR')
 GO TO 9999
 999 CALL mesage (-3,qhhl,NAME)
 9999 CALL mesage (-61,0,NAME)
 1000 CALL CLOSE (qhhl,1)
 RETURN
END SUBROUTINE frd2i
