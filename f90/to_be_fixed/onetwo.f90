SUBROUTINE onetwo(*,ix,x,dx,itermm)
!*******
!     PROGRAM TO SOLVE A MATRIX OF ORDER ONE OR TWO FOR DECOMP
!*******
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: ix(1)
 REAL, INTENT(IN OUT)                     :: x(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dx(6)
 INTEGER, INTENT(IN OUT)                  :: itermm
 DOUBLE PRECISION :: det,mindia,dz  ,da
 INTEGER :: sysbuf,rdp,dum
 INTEGER :: typel
 INTEGER :: scrflg,jposl,bbar,cbcnt,r,bbbar1 ,bbbar,sr2fl ,sr2fil
 INTEGER :: rd,wrt,rew,eofnrw
 DIMENSION sub(2)
 
 COMMON /system/sysbuf
 COMMON /dcompx/ifila(7),ifill(7),ifilu(7),dum(3),det,power, nx,mindia
 COMMON /names/ rd,rdrew,wrt,wrtrew,rew,norew,eofnrw ,rsp,rdp
 COMMON /zblpkx/dz(2),jj
 COMMON /packx/itype1,itype2,iy,jy,incry
 COMMON /unpakx/itypex,ixy,jxy,incrx
 
 EQUIVALENCE (ifila(2),ncol),(ifill(5),typel),(sr2fil,dum(2))
 
 DATA sub/4HONET,4HWO  /
 
! ----------------------------------------------------------------------
 
 ibuf1 = nx-sysbuf
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2-sysbuf
 ifile = ifilu(1)
 CALL CLOSE(dum(2),rew)
 IF(itermm == 1)ifile = dum(2)
 CALL gopen(ifile,ix(ibuf3),1)
 CALL gopen(ifila,ix(ibuf1),0)
 itypex = rdp
 itype1 = rdp
 itype2 = typel
 incrx = 1
 incry = 1
 IF(ncol == 2)GO TO 100
 IF( ncol /= 1)GO TO 5000
!*******
!     SOLVE A (1X1)
!*******
 ixy = 1
 jxy = 1
 CALL unpack(*5060,ifila(1),dx)
 det = dx(1)
 mindia = DABS(dx(1))
 iy = 1
 jy = 1
 CALL pack(dx,ifile,ifilu)
 dx(1) = 0.d0
 CALL pack(dx,ifill(1),ifill)
 IF(itermm == 0)GO TO 90
 CALL CLOSE(ifile,eofnrw)
 GO TO 95
 90 CALL CLOSE(ifile,rew)
 95 CALL CLOSE(ifila(1),rew)
 CALL CLOSE(ifill(1),rew)
 RETURN
 100 ixy = 1
!*******
!     SOLVE A (2X2)
!*******
 jxy = 2
 CALL unpack(*5060,ifila(1),dx)
 CALL unpack(*5060,ifila(1),dx(3))
 a = 1.
 IF(DABS(dx(1)) >= DABS(dx(2)))GO TO 150
!*******
!     PERFORM INTERCHANGE
!*******
 det = dx(1)
 dx(1) = dx(2)
 dx(2) = det
 det = dx(3)
 dx(3) = dx(4)
 dx(4) = det
 a = -1.
 150 CONTINUE
 dx(2) = dx(2)/dx(1)
 dx(4) = dx(4)-dx(2)*dx(3)
 det = dx(4)*dx(1)*a
 IF(dx(1) == 0.d0 .OR. dx(4) == 0.d0)GO TO 5060
 mindia = DMIN1 (DABS(dx(1)),DABS(dx(4)))
 iy = 1
 jy = 2
 dx(5) = 0.0D0
 IF(a < 0.0)  dx(5) = 1.0D0
 dx(6) = dx(2)
 CALL pack(dx(5),ifill(1),ifill)
 dx(6) = 0.
 jy = 1
 CALL pack(dx(6),ifill(1),ifill)
 IF(itermm == 1)GO TO 160
 dx(2) = dx(3)
 dx(3) = dx(4)
 dx(4) = dx(2)
 jy = 2
 CALL pack(dx(3),ifile,ifilu)
 iy = 2
 CALL pack(dx,ifile,ifilu)
 GO TO 90
 160 jy = 1
 CALL pack(dx,ifile,ifilu)
 jy=2
 CALL pack(dx(3),ifile,ifilu)
 CALL CLOSE(ifile,eofnrw)
 GO TO 95
 ENTRY finwrt(iterm,scrflg,sr2fl,jposl,i1sp,bbar,i1,cbcnt,  &
     ipak,r,bbbar1,bbbar,i6sp,i4,i4sp,ix,dx,x,lcol)
 ibuf1 = nx-sysbuf
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2-sysbuf
 CALL CLOSE(ifila(1),rew)
 CALL gopen(sr2fil,ix(ibuf1),wrt)
 CALL CLOSE(sr2fil,eofnrw)
 k=0
 CALL gopen(ifill,ix(ibuf2),wrt)
 IF(scrflg == 0)GO TO 2005
 CALL gopen(sr2fl,ix(ibuf3),rd)
 2005 ll = 0
 2010 jposl = jposl+1
 CALL bldpk(rdp,typel,ifill(1),0,0)
 in1 = i1sp+k
 jj = jposl
 dz(1) = ix(in1)
 CALL zblpki
 kk = 0
 iend = MIN0(bbar,ncol-jj)
 IF(iend == 0)GO TO 2030
 in1 = i1+ll*bbar
 2020 jj = jj+1
 in2 = in1+kk
 dz(1) =dx(in2)
 CALL zblpki
 kk = kk+1
 IF(kk-iend < 0) THEN
   GO TO  2020
 ELSE IF (kk-iend == 0) THEN
   GO TO  2030
 ELSE
   GO TO  5050
 END IF
 2030 IF(cbcnt == 0)GO TO 2050
!*******
!     PACK ACTIVE ROW ELEMENTS ALSO
!*******
 kk = 0
 2035 in1 = i6sp + kk
 in2 = i4 + ix(in1)*bbbar + k
 dz(1) = dx(in2)
 IF(dz(1) == 0.d0)GO TO 2040
 in1 = i4sp + ix(in1)
 jj = ix(in1)
 CALL zblpki
 2040 kk = kk + 1
 IF(kk < cbcnt)GO TO 2035
 2050 CALL bldpkn(ifill(1),0,ifill)
 ll = ll + 1
 k = k + 1
 IF(k == lcol)GO TO 2080
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
 2070 ll =ll-1
 in1 = i1+ll*bbar
 CALL fread(sr2fl,dx(in1),2*bbar,0)
 GO TO 2010
 2080 CALL CLOSE(ifill(1),rew)
 IF(scrflg > 0)CALL CLOSE(sr2fl,rew)
 IF(iterm /= 0)RETURN
!*******
!     RE-WRITE THE UPPER TRIANGLE WITH THE RECORDS IN THE REVERSE ORDER
!*******
 incrx = 1
 incry = 1
 itype1 = typel
 itype2 = typel
 itypex = typel
 ifilu(2) = 0
 ifilu(6) = 0
 ifilu(7) = 0
 CALL gopen(sr2fil,ix(ibuf1),rd)
 CALL gopen(ifilu,ix(ibuf2),1)
 DO  i = 1,ncol
   ixy = 0
   CALL bckrec(sr2fil)
   CALL unpack(*5060,sr2fil,ix)
   CALL bckrec(sr2fil)
   kk = jxy-ixy+1
   k = kk/2
   kk = kk + 1
   IF(typel == 1)GO TO 2095
   DO  j = 1,k
     l = kk-j
     da = dx(j)
     dx(j) = dx(l)
     dx(l) = da
   END DO
   GO TO 2100
   2095 DO  j = 1,k
     l = kk-j
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
 5050 no = -25
 GO TO 5500
 5060 RETURN 1
 5500 CALL mesage(no,0,sub)
 RETURN
END SUBROUTINE onetwo
