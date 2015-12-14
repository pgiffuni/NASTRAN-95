SUBROUTINE gust3 (qhjk,wj,pp,gustl,pdel,pgust,q,nfreq,nload,  &
        nrowj,ncolw)
     
!     THE PURPOSE OF THIS ROUTINE IS TO MULTIPLY QHJK(+) BY WJ
!     FORMING PDEL
!     PDEL IS THEN MULTIPLIED BY  Q*WG*PP(W)  FORMING PGUST
 
 
 INTEGER, INTENT(IN)                      :: qhjk
 INTEGER, INTENT(IN OUT)                  :: wj
 INTEGER, INTENT(IN OUT)                  :: pp
 INTEGER, INTENT(IN OUT)                  :: gustl
 INTEGER, INTENT(IN OUT)                  :: pdel
 INTEGER, INTENT(IN OUT)                  :: pgust
 REAL, INTENT(IN)                         :: q
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER, INTENT(IN)                      :: nload
 INTEGER, INTENT(IN)                      :: nrowj
 INTEGER, INTENT(IN OUT)                  :: ncolw
 INTEGER :: iz(1),sysbuf,mcb(7), NAME(2)
 COMMON /packx / itc1,itc2,ii1,jj1,incr1
 COMMON /unpakx/ itc,ii,jj,incr
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (iz(1),z(1))
 DATA    NAME  / 4HGUST,1H3 /
 
!     INITIALIZE
 
 ibuf1 = korsz(iz)- sysbuf+1
 ibuf2 = ibuf1- sysbuf
 ibuf3 = ibuf2- sysbuf
 incr1 = 1
 incr  = 1
 ibuf4 = ibuf3- sysbuf
 mcb(1)= qhjk
 CALL rdtrl(mcb)
 itc  = 3
 itc1 = itc
 itc2 = itc
 CALL gopen (wj,iz(ibuf1),0)
 CALL gopen (qhjk,iz(ibuf2),0)
 CALL gopen (pdel,iz(ibuf3),1)
 
!     SET UP TO PACK
 
 it1 = 1
 jj1 = mcb(3) / nrowj
 nrqhj = mcb(3)
 ntqhj = nrqhj*2
 CALL makmcb (mcb,pdel,jj1,2,itc2)
 ii   = 1
 iqhj = 2*nfreq+1
 iwj  = iqhj+ntqhj
 ntwz = nrowj*2
 ipdel= iwj + ntwz
 ntpdel = jj1*2
 nz = ibuf4-1 - ipdel + 2*jj1
 IF (nz  < 0) CALL mesage (-8,0,NAME)
 DO  i = 1,nfreq
   jj = nrqhj
   CALL unpack (*10,qhjk,z(iqhj))
   
!     MULTIPY EACH IMAGINARY PART BY K
   
   DO  j = 1,ntqhj,2
     z(iqhj+j) = z(iqhj+j)*z(2*i)
   END DO
   GO TO 20
   
!     NULL COLUMN
   
   10 CALL zeroc (z(iqhj),ntqhj)
   20 CONTINUE
   
!     BRING WJ COLUMN INTO CORE
   
   jj = nrowj
   CALL unpack (*30,wj,z(iwj))
   GO TO 40
   30 CALL zeroc (z(iwj),ntwz)
   40 CONTINUE
   
!     MULTIPLY
   
   CALL gmmatc (z(iqhj),jj1,nrowj,0,z(iwj),nrowj,1,0,z(ipdel))
   CALL pack (z(ipdel),pdel,mcb)
 END DO
 CALL CLOSE (wj,1)
 CALL CLOSE (qhjk,1)
 CALL CLOSE (pdel,1)
 CALL wrttrl (mcb)
 CALL dmpfil (-pdel,z,nz)
 
!     REPEATEDLY READ PDEL MULTIPLYING BY Q,WG, AND PP
 
 CALL gopen (pdel,iz(ibuf1),0)
 CALL gopen (pp,iz(ibuf2),0)
 CALL gopen (gustl,iz(ibuf3),0)
 CALL gopen (pgust,iz(ibuf4),1)
 CALL makmcb (mcb,pgust,mcb(3),mcb(4),mcb(5))
 DO  i = 1,nload
   CALL REWIND (pdel)
   CALL skprec (pdel,1)
   CALL fread (gustl,iz,5,1)
   iz2 = 2
   qwg = q*z(iz2+1)
   DO  j = 1,nfreq
     jj = 1
     CALL unpack (*310,pp,z)
     qwgr = qwg * z(1)
     qwgc = qwg * z(iz2)
     GO TO 320
     310 CONTINUE
     qwgr = 0.0
     qwgc = 0.0
     320 CONTINUE
     jj = jj1
     CALL unpack (*330,pdel,z)
     GO TO 340
     330 CALL zeroc (z,ntpdel)
     340 CONTINUE
     DO  m = 1,ntpdel,2
       pgr = qwgr*z(m  ) - qwgc*z(m+1)
       pgc = qwgr*z(m+1) + qwgc*z(m  )
       z(m  ) = pgr
       z(m+1) = pgc
     END DO
     CALL pack (z,pgust,mcb)
   END DO
 END DO
 CALL CLOSE (pdel,1)
 CALL CLOSE (pp,1)
 CALL CLOSE (gustl,1)
 CALL CLOSE (pgust,1)
 CALL wrttrl (mcb)
 RETURN
END SUBROUTINE gust3
