SUBROUTINE com12  (*,ix,x,dx,itermm)
     
!*******
!      PROGRAM TO SOLVE A MATRIX OF ORDER ONE OR TWO FOR CDCOMP
!*******
 
 , INTENT(IN)                             :: *
 INTEGER, INTENT(IN)                      :: ix(1)
 REAL, INTENT(IN OUT)                     :: x(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dx(12)
 INTEGER, INTENT(IN OUT)                  :: itermm
 DOUBLE PRECISION :: det,mindia,dz,da
 INTEGER :: sysbuf,rdp,dum
 INTEGER :: typel
 INTEGER :: cdp
 INTEGER :: scrflg,jposl,bbar,cbcnt,r,bbbar1,bbbar, sr2fl,sr2fil
 DIMENSION       sub(2)
 COMMON /system/ sysbuf
 COMMON /cdcmpx/ ifila(7),ifill(7),ifilu(7),dum(3),det(2),power, nx,mindia
 COMMON /names / rd,rdrew,wrt,wrtrew,rew,norew,eofnrw ,rsp,rdp, csp,cdp
 COMMON /zblpkx/ dz(2),jj
 COMMON /packx / itype1,itype2,iy,jy,incry
 COMMON /unpakx/ itypex,ixy,jxy,incrx
 EQUIVALENCE     (ifila(2),ncol),(ifill(5),typel),(sr2fil,dum(2))
 DATA    sub(1), sub(2)  / 4HCOM1,4H2    /
 
 ibuf1 = nx-sysbuf
 ibuf2 = ibuf1-sysbuf
 CALL CLOSE(sr2fil,rew)
 ibuf3 = ibuf2-sysbuf
 ifile = ifilu(1)
 IF(itermm == 1) ifile = sr2fil
 CALL gopen(ifile,ix(ibuf3),wrtrew)
 CALL gopen(ifila(1),ix(ibuf1),rdrew)
 itypex = cdp
 itype1 = cdp
 itype2 = typel
 incrx  = 1
 incry  = 1
 IF(ncol == 2) GO TO 100
 IF(ncol /= 1) GO TO 5000
!*******
!     SOLVE A (1X1)
!*******
 ixy = 1
 jxy = 1
 CALL unpack(*5060,ifila(1),dx)
 det(1) = dx(1)
 det(2) = dx(2)
 mindia = DSQRT(dx(1)**2+dx(2)**2)
 iy = 1
 jy = 1
 CALL pack (dx,ifile,ifilu)
 dx(1) = 0.d0
 dx(2) = 0.d0
 CALL pack (dx,ifill(1),ifill)
 90 CALL CLOSE(ifile,rew)
 95 CALL CLOSE(ifila(1),rew)
 CALL CLOSE(ifill(1),rew)
 RETURN
 100 ixy = 1
!*******
!     SOLVE A (2X2)
!*******
 jxy = 2
 CALL unpack(*5060,ifila(1),dx   )
 CALL unpack(*5060,ifila(1),dx(5))
 a = 1.
 IF(dx(1)**2+dx(2)**2 >= dx(3)**2+dx(4)**2) GO TO 150
!*******
!     PERFORM INTERCHANGE
!*******
 det(1) = dx(1)
 dx(1)  = dx(3)
 dx(3)  = det(1)
 det(1) = dx(2)
 dx(2)  = dx(4)
 dx(4)  = det(1)
 det(1) = dx(5)
 dx(5)  = dx(7)
 dx(7)  = det(1)
 det(1) = dx(6)
 dx(6)  = dx(8)
 dx(8)  = det(1)
 a      = -1.
 150 CONTINUE
 det(1) = (dx(3)*dx(1)+dx(4)*dx(2))/(dx(1)**2+dx(2)**2)
 dx(4)  = (dx(4)*dx(1)-dx(3)*dx(2))/(dx(1)**2+dx(2)**2)
 dx(3)  = det(1)
 dx(7)  = dx(7)-dx(3)*dx(5)+dx(4)*dx(6)
 dx(8)  = dx(8)-dx(3)*dx(6)-dx(4)*dx(5)
 det(1) = dx(1)*dx(7)-dx(2)*dx(8)
 det(2) = dx(2)*dx(7)+dx(1)*dx(8)
 IF((dx(1) == 0.d0 .AND. dx(2) == 0.d0) .OR. (dx(7) == 0.d0.AND.  &
     dx(8) == 0.d0)) GO TO 5060
 mindia = DMIN1(DSQRT(dx(1)**2+dx(2)**2),DSQRT(dx(7)**2+dx(8)**2))
 iy = 1
 jy = 2
 dx( 9) = 0.0D0
 dx(10) = 0.0D0
 IF(a < 0.) dx(9) = 1.d0
 dx(11) = dx(3)
 dx(12) = dx(4)
 CALL pack(dx(9),ifill(1),ifill)
 dx( 9) = 0.d0
 dx(10) = 0.d0
 jy = 1
 CALL pack(dx(9),ifill(1),ifill)
 IF(itermm == 1) GO TO 160
 dx(3) = dx(5)
 dx(5) = dx(7)
 dx(7) = dx(3)
 dx(3) = dx(6)
 dx(6) = dx(8)
 dx(8) = dx(3)
 jy    = 2
 CALL pack(dx(5),ifilu(1),ifilu)
 iy = 2
 CALL pack (dx,ifilu(1),ifilu)
 GO TO 90
 160 CALL pack(dx,ifile,ifilu)
 jy = 2
 CALL pack(dx(5),ifile,ifilu)
 CALL CLOSE(ifile,eofnrw)
 GO TO 95
 
 
 ENTRY comfin (iterm,scrflg,sr2fl,jposl,i1sp,bbar,i1,cbcnt,  &
     ipak,r,bbbar1,bbbar,i6sp,i4,i4sp,ix,dx,x,lcol)
 
 ibuf1 = nx-sysbuf
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2-sysbuf
 CALL CLOSE(ifila(1),rew)
 CALL OPEN(*5010,sr2fil,ix(ibuf1),wrt)
 CALL CLOSE(sr2fil,eofnrw)
 k=0
 NAME =  ifill(1)
 CALL OPEN(*5010,ifill(1),ix(ibuf2),wrt)
 IF(scrflg == 0) GO TO 2005
 NAME = sr2fl
 CALL OPEN(*5010,sr2fl,ix(ibuf3),rd)
 2005 ll = 0
 2010 jposl = jposl+1
 CALL bldpk(cdp,typel,ifill(1),0,0)
 in1 = i1sp+k
 jj  = jposl
 dz(1) = ix(in1)
 dz(2) = 0.d0
 CALL zblpki
 kk   = 0
 iend = MIN0(bbar,ncol-jj)
 IF(iend == 0) GO TO 2030
 in1 = i1 +ll*bbar*2
 2020 jj  = jj+1
 in2 = in1+kk+kk
 dz(1) = dx(in2)
 dz(2) = dx(in2+1)
 CALL zblpki
 kk = kk+1
 IF(kk-iend < 0) THEN
   GO TO  2020
 ELSE IF (kk-iend == 0) THEN
   GO TO  2030
 ELSE
   GO TO  5050
 END IF
 2030 IF(cbcnt == 0) GO TO 2050
!*******
!     PACK ACTIVE ROW ELEMENTS ALSO
!*******
 kk  = 0
 2035 in1 = i6sp + kk
 in2 = i4 +(ix(in1)*bbbar+k)*2
 dz(1) = dx(in2)
 dz(2) = dx(in2+1)
 IF(dz(1) == 0.d0 .AND. dz(2) == 0.d0) GO TO 2040
 in1 = i4sp + ix(in1)
 jj  = ix(in1)
 CALL zblpki
 2040 kk = kk + 1
 IF(kk < cbcnt) GO TO 2035
 2050 CALL bldpkn(ifill(1),0,ifill)
 ll = ll + 1
 k  = k + 1
 IF (k == lcol) GO TO 2080
 IF(k-r+1 < 0) THEN
   GO TO  2010
 ELSE IF (k-r+1 == 0) THEN
   GO TO  2060
 ELSE
   GO TO  2070
 END IF
 2060 IF(r-bbbar1 < 0.0) THEN
   GO TO  2070
 ELSE IF (r-bbbar1 == 0.0) THEN
   GO TO  2010
 ELSE
   GO TO  5050
 END IF
 2070 ll  = ll-1
 in1 = i1+ll*bbar*2
 ibbar4 = 4 * bbar
 CALL READ(*5020,*5030,sr2fl,dx(in1),ibbar4,0,no)
 GO TO 2010
 2080 CALL CLOSE(ifill(1),rew)
 IF(scrflg > 0)CALL CLOSE(sr2fl,rew)
 IF(iterm /= 0) RETURN
!*******
!     RE-WRITE THE UPPER TRIANGLE WITH THE RECORDS IN THE REVERSE ORDER
!*******
 incrx  = 1
 incry  = 1
 itype1 = typel
 itype2 = typel
 itypex = typel
 ifilu(2) = 0
 ifilu(6) = 0
 ifilu(7) = 0
 NAME = sr2fil
 CALL OPEN(*5010,sr2fil,ix(ibuf1),rd)
 CALL gopen(ifilu(1),ix(ibuf2),wrtrew)
 DO  i = 1,ncol
   ixy = 0
   CALL bckrec(sr2fil)
   CALL unpack(*5060,sr2fil,ix)
   CALL bckrec(sr2fil)
   kk = jxy-ixy+1
   k  = kk/2
   kk = kk + 1
   IF(typel == 1) GO TO 2095
   IF(typel == 4) GO TO 2061
   DO  j = 1,k
     l  = kk-j
     da = dx(j)
     dx(j) = dx(l)
     dx(l) = da
   END DO
   GO TO 2100
   2061 kk = kk+kk-1
   k  = k+k
   DO  j = 1,k,2
     l  = kk-j-1
     da = dx(l)
     dx(l) = dx(j)
     dx(j) = da
     da = dx(l+1)
     dx(l+1) = dx(j+1)
     dx(j+1) = da
   END DO
   GO TO 2100
   2095 DO  j = 1,k
     l    = kk-j
     a    = x(j)
     x(j) = x(l)
     x(l) = a
   END DO
   2100 iy = ncol-jxy+1
   jy = ncol-ixy+1
   CALL pack(ix,ifilu(1),ifilu)
 END DO
 CALL CLOSE(ifilu(1),rew)
 CALL CLOSE(sr2fil,rew)
 RETURN
 5000 no = -8
 GO TO 5500
 5010 no = -1
 GO TO 5500
 5020 no = -2
 GO TO 5500
 5030 no = -3
 GO TO 5500
 5050 no = -25
 GO TO 5500
 5060 RETURN 1
 5500 CALL mesage(no,NAME,sub(1))
 RETURN
END SUBROUTINE com12
