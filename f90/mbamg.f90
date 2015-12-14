SUBROUTINE mbamg (INPUT,ajjl,skj)
     
!     DRIVER FOR MACH BOX THEORY
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: ajjl
 INTEGER, INTENT(IN OUT)                  :: skj
 LOGICAL :: cntrl2,cntrl1,crank1,crank2,asym
 INTEGER :: sysbuf, NAME(2),iz(1),buf1,scr2
 REAL :: mach
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /mboxa / x(12),y(12),tang(10),ang(10),cotang(10)
 COMMON /mboxc / njj ,crank1,crank2,cntrl1,cntrl2,nbox,  &
     npts0,npts1,npts2,asym,gc,cr,mach,beta,ek,ekbar,  &
     ekm,boxl,boxw,boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2, kct,kc1t,kc2t
 COMMON /system/ sysbuf,nout
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /packx / iti,it0,ii,nn,incr
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    NAME  / 4HMBAM,4HG    /
 DATA    nhcore, nhcapf,nhcont /4HCORE,4HCAPF,4HCONT/
 DATA    scr2  / 302 /
 
!     SCR2 CONTAINS THE INTERPOLATED POINTS
 
!     2 * KCT FOR NPTS0 POINTS
!     2 * KC1T FOR NPTS1 POINTS
!     2 * KC2T FOR NPTS2 POINTS
 
 
!     OPEN CORE POINTERS FIXED DIMENSIONS
 
 nw1   = 1
 nwn   = 51
 nc21  = 101
 nc2n  = 151
 nc1   = 201
 ncn   = 251
 nd1   = 301
 ndn   = 351
 nxk   = 401
 nyk   = 601
 nxk1  = 801
 nyk1  = 926
 nxk2  = 1051
 nyk2  = 1176
 nxwte = 1301
 nywte = 1351
 nkte  = 1401
 nkte1 = 1451
 nkte2 = 1501
 nparea=1551
 icorr = 9051
 
!     INITITALIZE  PUT HEADER DATA IN MBOXC
 
 icore = korsz(iz) - 4*sysbuf
 buf1  = icore - sysbuf
 CALL fread (INPUT,njj,9,0)
 asym = .false.
 IF( nd == -1 ) asym = .true.
 mach = fmach
 beta = SQRT((mach*mach)-1.0)
 CALL fread (INPUT,z,24,0)
 
!     MOVE X AND Y TO MBOXA
 
 l = 0
 DO  i = 1,23,2
   l = l + 1
   x(l)  = z(i)
   y(l)  = z(i+1)
 END DO
 CALL mbgeod
 ek    = (2.0*cr/refc)*rfk
 cmax  =  AMAX1(x(4),x(5),x(6))
 boxl  =  cmax/(FLOAT(nbox) + 0.50)
 boxw  =  boxl/beta
 nsb   =  y(3)/boxw + 0.5
 nsb   =  MIN0(nsb,50)
 boxw  =  y(3)/(FLOAT(nsb) - 0.50)
 boxl  =  boxw*beta
 ncb   =  cmax/boxl + 0.999
 
!     CALL MBREG TO GENERATE BOXES
 
 icrq = icorr - buf1
 IF (icorr > buf1) GO TO 996
 20 CALL mbreg (ireg,z(nw1),z(nwn),z(nc21),z(nc2n),z(nc1),z(ncn),  &
     z(nd1),z(ndn),z(nxk),z(nyk),z(nxk1),z(nyk1),z(nxk2),  &
     z(nyk2),z(nxwte),z(nywte),z(nkte),z(nkte1),z(nkte2), z(nparea))
 IF (ireg /= 2) GO TO 30
 IF (nbox < 2) GO TO 999
 nbox = nbox - 1
 GO TO 20
 30 CALL mbplot (z(nw1),z(nd1),z(nwn),z(nc21),z(nc2n),z(nc1), z(ncn),z(ndn))
 
!     CALL MBMODE TO GENERATE MODE LIKE DATA
 
 CALL gopen (scr2,z(buf1),1)
 CALL mbmode (INPUT,scr2,icorr,buf1,z,npts0,kct,z(nxk),z(nyk),is, cr)
 IF (is == 2) GO TO 997
 IF (cntrl1) CALL mbmode (INPUT,scr2,icorr,buf1,z,npts1,kc1t,  &
     z(nxk1),z(nyk1),is,cr)
 IF (is == 2) GO TO 997
 IF (cntrl2) CALL mbmode (INPUT,scr2,icorr,buf1,z,npts2,kc2t,  &
     z(nxk2),z(nyk2),is,cr)
 IF (is == 2) GO TO 997
 CALL CLOSE (scr2,1)
 ekbar = (ek*boxl*mach*mach)/(beta*beta)
 ekm   = ekbar/mach
 CALL fread (INPUT,0,0,1)
 CALL bug (nhcore,80,z,nyk1-1)
 CALL bug (nhcore,80,z(nyk1),nparea-nyk1)
 CALL dmpfil (scr2 ,z(icorr),buf1-icorr)
 
!     MORE DIMENSIONS
 
 IF (MOD(icorr,2) == 0) icorr = icorr + 1
 ncap  = icorr
 ncaph = ncb*(ncb+1)/2
 
!     COMPLEX PHIS
 
 icorr = ncap + ncaph*2
 icrq  = icorr - buf1
 IF (icorr > buf1) GO TO 996
 CALL mbcap (ncaph,z(ncap))
 icorr = ncap + ncaph*2
 CALL bug (nhcapf,80,z(ncap),ncaph*2)
 
!     PUT OUT SKJ
 
 iti = 1
 it0 = 3
 ii  = isk
 nsk = nsk + 1
 nn  = nsk
 rm  = 1.0
 DO  i = 1,njj
   CALL pack (rm,skj,tskj)
   ii  = ii + 1
   IF (i == njj) CYCLE
   nn  = nn + 1
 END DO
 isk = ii
 nsk = nn
 
!     SET UP FOR COLUMN OF AJJL
 
 iti = 3
 it0 = 3
 ii  = nrow + 1
 nn  = nrow + njj
 
!     GET AJJL MATRIX TERMS
!     MORE DIMENSIONS
 
 nphit = icorr
 ndss  = nphit + (3*nsbd)*2
 nq    = ndss  + (ncb*nsbd)*2
 nq1   = nq + kct*2
 nq2   = nq1 + kc1t*2
 na    = nq2 + kc2t*2
 icorr = na + njj*2
 CALL bug (nhxect,100,x,54)
 CALL bug (nhcont,100,njj,30)
 icrq  = icorr - buf1
 IF (icorr > buf1) GO TO 996
 CALL mbdpdh (ajjl,z(nxk),z(nyk),z(nxk1),z(nyk1),z(nxk2),z(nyk2),  &
     z(nxwte),z(nywte),z(nparea),z(ncap),z(nphit),z(ndss),  &
     z(nq),z(nq1),z(nq2),z(ndn),z(nd1),z(nw1),z(nwn),  &
     z(nkte),z(nkte1),z(nkte2),z(nc1),ncb,nsbd,scr2, z(buf1),z(na))
 nrow = nrow + njj
 1000 RETURN
 
!     ERROR MESSAGES
 
 997 WRITE  (nout,9971) ufm
 9971 FORMAT (a23,' 2424, MACH BOX CONTROL POINTS IMPROPER SINGULAR ',  &
     'MATRIX RESULTED')
 GO TO 998
 999 WRITE  (nout,9991) ufm
 9991 FORMAT (a23,' 2425, MACH BOX GENERATION OF BOXES FAILED')
 998 CALL mesage (-37,0,NAME)
 996 CALL mesage (-8,icrq,NAME)
 GO TO 1000
END SUBROUTINE mbamg
