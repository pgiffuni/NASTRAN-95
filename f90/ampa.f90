SUBROUTINE ampa (aero,qjh,qhh,ajjl,qhhlo,qjhlo,INDEX,imax,iany)
     
!     THE PURPOSE OF THIS ROUTINE IS TO
!         1. INITIALIZE QHJ AND QHH
!         2. COPY USEFUL DATA FROM QJH AND QHH TO QHHLO AND  QJHLO
!         3. SET UP INDEX, IMAX,IANY, AND AMPCOM
 
!     OPEN CORE IS LAID OUT AS FOLLOWS
 
!     CONTENTS                POINTER                  LENGTH
!     --------                -------                  ------
!     AJJL HEADER
!     NCOL
!     NSUB
!     M-K PAIRS               IAJJL                    2*NSUB +2
!      .
!      .
!      .
!     AERO RECORD 2           IAERO                    2* IMAX
!     M- K  PAIRS
!      .
!      .
!      .
!     QHH HEADER RECORD(RST)
!     NOH (OLD)
!     M- K PAIRS              IQHH                     2*NQHH
!       .
!       .
!       .
!     BUFFER2                 IBUF2
!     BUFFER1                 IBUF1
 
 
!     SPECIAL CODE EXISTS IN CASE AJJK HEADER HAS ONLY 2 WORDS
 
 
 INTEGER, INTENT(IN)                      :: aero
 INTEGER, INTENT(IN)                      :: qjh
 INTEGER, INTENT(IN)                      :: qhh
 INTEGER, INTENT(IN)                      :: ajjl
 INTEGER, INTENT(IN)                      :: qhhlo
 INTEGER, INTENT(IN)                      :: qjhlo
 INTEGER, INTENT(IN OUT)                  :: INDEX
 INTEGER, INTENT(OUT)                     :: imax
 INTEGER, INTENT(OUT)                     :: iany
 INTEGER :: xqhhl,sysbuf, NAME(2),mcbajj(7),FILE,ajjcol,qhhcol
 REAL :: z(1)
 COMMON /unpakx/ it1,ii,jj,incr
 COMMON /packx / it2,it3,ii1,jj1,incr1
 COMMON /system/ sysbuf,nout,skp(52),iprec
 COMMON /BLANK / noue,xqhhl,igust
 COMMON /ampcom/ ncol,nsub,xm,xk,ajjcol,qhhcol,ngp,ngpd(2,30),  &
     mcbqhh(7),mcbqjh(7),noh,idjh,mcbrjh(7)
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    NAME  / 4HAMPA,4H..../
 
!     INITIALIZE
 
 mcbajj(1) = ajjl
 CALL rdtrl (mcbajj)
 mcbqhh(1) = qhh
 CALL rdtrl (mcbqhh)
 mcbqjh(1) = qjh
 CALL rdtrl (mcbqjh)
 iany = 1
 ibuf1= korsz(iz) - sysbuf + 1
 ibuf2= ibuf1 - sysbuf
 
!     EXTRACT DATA FROM AJJL HEADER
 
 FILE = ajjl
 CALL OPEN (*900,ajjl,iz(ibuf1),0)
 CALL READ (*910,*920,ajjl,iz,-2,0,iflag)
 CALL READ (*910,*30,ajjl,iz,ibuf2-1,0,iflag)
 GO TO 980
 
!     PROCESS AJJL  DATA
 
 30 CALL CLOSE (ajjl,1)
 noajjh= 0
 iajjl = 4
 IF (iflag == 0) GO TO 50
 ncol = iz(1)
 izx  = 3
 nsub = MIN0(iz(izx),mcbajj(2)/ncol)
 ngp  = iz(2*nsub+4)
 k    = 2*nsub + 5
 iaero= k - 1
 DO  i = 1,ngp
   ngpd(1,i) = iz(k  )
   ngpd(2,i) = iz(k+1)
   k = k + 3
 END DO
 GO TO 55
 
!     NO AJJ HEADER DATA
 
 50 noajjh = 1
 iaero  = 3
 
!     BRING IN AERO DATA
 
 55 CONTINUE
 CALL gopen (aero,iz(ibuf1),0)
 FILE = aero
 CALL fwdrec (*910,aero)
 nz = ibuf2 - iaero
 CALL READ (*910,*60,aero,iz(iaero),nz,0,iflag)
 GO TO 980
 
!     AERO DATA IN CORE
 
 60 CALL CLOSE (aero,1)
 imax = iflag/2
 IF (noajjh == 0) GO TO 70
 
!     FIX UP FOR AJJ MISSING HEADER
 
 ncol = mcbajj(2)/imax
 nsub = imax
 ngp  = 1
 ngpd(1,1) = 1
 ngpd(2,1) = ncol
 iaero = iflag + 3
 k = iaero
 DO  i = 1, iflag
   iz(k) = iz(i+2)
   k = k+1
 END DO
 
!     PUT HEADERS FROM OLD QHH IN CORE
 
 70 IF (xqhhl == 1) GO TO 80
 FILE = qhh
 CALL OPEN (*900,qhh,iz(ibuf1),0)
 CALL fread (qhh,iz,-2,0)
 iqhh = iaero + 2*imax + 2
 nz   = nz - 2*imax
 CALL READ (*910,*75,qhh,iz(iqhh),nz,0,iflag)
 GO TO 980
 75 CALL CLOSE (qhh,1)
 iqhh = iqhh + 2
 nqhh = MIN0((iflag-2)/2,mcbqhh(2)/noh)
 
!     BUILD INDEX FILE
 
 80 CONTINUE
 i = 0
 CALL gopen (INDEX,iz(ibuf1),1)
 90 CONTINUE
 xm =  z(iaero+i  )
 xk =  z(iaero+i+1)
 
!     SEARCH FOR COLUMN NUMBER IN AJJL
 
 j  =  0
 100 CONTINUE
 xma = z(iajjl+j  )
 xka = z(iajjl+j+1)
 IF (xma == xm .AND. xka == xk) GO TO 120
 j =  j + 2
 IF (j >= 2*nsub) CALL mesage (-7,0,NAME)
 GO TO 100
 
!     FOUND IN AJJL
 
 120 CONTINUE
 ajjcol = (j/2)*ncol + 1
 
!     SEARCH FOR COLUMN NUMBER IN QHH
 
 qhhcol = 0
 IF (xqhhl == 1) GO TO 140
 j = 0
 130 CONTINUE
 xma = z(iqhh+j  )
 xka = z(iqhh+j+1)
 IF (xma == xm .AND. xka == xk) GO TO 150
 j = j + 2
 IF (j >= 2*nqhh) GO TO 140
 GO TO 130
 
!     FOUND IN QHH
 
 150 qhhcol = (j/2)*noh + 1
 
!     WRITE ON INDEX
 
 140 CALL WRITE (INDEX,xm,4,1)
 IF (qhhcol == 0) iany = 0
 i = i + 2
 IF (i >= 2*imax) GO TO 200
 GO TO 90
 
!     DONE WITH INDEX
 
 200 CALL CLOSE (INDEX,1)
 
!     COPY OLD  QHH  ONTO QHHLO
 
 IF (xqhhl == 1) GO TO 300
 it1  = mcbqhh(5)
 it2  = it1
 it3  = it1
 incr = 1
 incr1= 1
 IF (mcbqhh(1) <= 0) GO TO 230
 CALL gopen (qhh,iz(ibuf1),0)
 CALL gopen (qhhlo,iz(ibuf2),1)
 nclqhh    = mcbqhh(2)
 mcbqhh(2) = 0
 mcbqhh(6) = 0
 mcbqhh(7) = 0
 mcbqhh(1) = qhhlo
 CALL cyct2b (qhh,qhhlo,nclqhh,iz,mcbqhh)
 CALL CLOSE (qhh,1)
 CALL CLOSE (qhhlo,1)
 CALL wrttrl (mcbqhh)
 
!     COPY OLD QJH ONTO QJHLO
 
 230 CONTINUE
 
!     COPY QJH ONTO QJHLO
 
 IF (mcbqjh(1) <= 0) GO TO 250
 CALL gopen (qjh,iz(ibuf1),0)
 CALL gopen (qjhlo,iz(ibuf2),1)
 nclqjh    = mcbqjh(2)
 mcbqjh(1) = qjhlo
 mcbqjh(2) = 0
 mcbqjh(6) = 0
 mcbqjh(7) = 0
 CALL cyct2b (qjh,qjhlo,nclqjh,iz,mcbqjh)
 CALL CLOSE (qjh,1)
 CALL CLOSE (qjhlo,1)
 CALL wrttrl (mcbqjh)
 250 CONTINUE
 
!     PUT HEADERS ON NEW OUTPUT FILES
 
 300 CONTINUE
 IF (mcbqhh(1) <= 0) GO TO 350
 FILE  = qhh
 CALL OPEN (*900,qhh,iz(ibuf1),1)
 CALL fname (qhh,mcbqhh)
 CALL WRITE (qhh,mcbqhh,2,0)
 CALL WRITE (qhh,noh,1,0)
 CALL WRITE (qhh,imax,1,0)
 CALL WRITE (qhh,iz(iaero),2*imax,1)
 CALL CLOSE (qhh,3)
 mcbqhh(1) = qhh
 mcbqhh(2) = 0
 mcbqhh(3) = noh
 mcbqhh(4) = 2
 mcbqhh(5) = 2 + iprec
 mcbqhh(6) = 0
 mcbqhh(7) = 0
 350 CONTINUE
 IF (mcbqjh(1) <= 0) GO TO 360
 FILE = qjh
 CALL OPEN (*900,qjh,iz(ibuf1),1)
 CALL fname (qjh,mcbqjh)
 CALL WRITE (qjh,mcbqjh,2,0)
 CALL WRITE (qjh,noh,1,0)
 CALL WRITE (qjh,imax,1,0)
 CALL WRITE (qjh,iz(iaero),2*imax,1)
 CALL CLOSE (qjh,3)
 mcbqjh(1) = qjh
 mcbqjh(2) = 0
 mcbqjh(3) = ncol
 mcbqjh(4) = 2
 mcbqjh(5) = 2 + iprec
 mcbqjh(6) = 0
 mcbqjh(7) = 0
 360 CONTINUE
 iany = 0
 
!     PUT HEADER ON QHJL
 
 IF (igust <= 0) RETURN
 FILE = mcbrjh(1)
 CALL OPEN (*900,FILE,iz(ibuf1),1)
 CALL fname (FILE,mcbrjh(2))
 CALL WRITE (FILE,mcbrjh(2),2,0)
 CALL WRITE (FILE,noh,1,0)
 CALL WRITE (FILE,imax,1,0)
 CALL WRITE (FILE,iz(iaero),2*imax,1)
 CALL CLOSE (FILE,3)
 CALL makmcb (mcbrjh,FILE,ncol,2,2+iprec)
 CALL wrttrl (mcbrjh)
 RETURN
 
!     ERROR MESSAGES
 
 900 ip1 = -1
 901 CALL mesage (ip1,FILE,NAME)
 910 ip1 = -2
 GO TO 901
 920 ip1 = -3
 GO TO 901
 980 ip1 = -8
 GO TO 901
END SUBROUTINE ampa
