SUBROUTINE em1d (eltype,istart,itype,ncount,ido,iwords,nbdys,  &
        all,nelout)
     
!     COMPUTE LOAD DUE TO MAGNETIC FIELD,  K*A + F = 0
!     SOLVE FOR -F
 
!     USE ELEMENT COORDINATES  FOR  ROD
 
!     SET UP COMMON BLOCKS, TABLES
 
!     KSYSTM(1) = 1ST POS. OF OPEN CORE
!     KSYSTM(2) = OUTPUT FILE NO.
!     KSYSTM(56) NE 0 FOR HEAT TRANSFER OPTION
 
!     Z      = OPEN CORE ARRAY
!     OUTPT  = OUTPUT FILE NO.
!     NELEMS = NO OF ELEMENTS (TYPES) IN THIS TABLE
!     LAST   = LOC OF 1ST WORD OF LAST ENTRY(EL) IN TABLE
!     INCR   = MAX NO WDS ALLOWED IN ANY ENTRY
 
!     BUF1   = BUFFER FOR EST
!     EST    = ELEMENT SUMMARY TABLE(PROG MAN 2.3.56)
!     SLT    = STAIC LOADS TABLE(2.3.51)
!     SYSTEM 2.4.13 PROG MANUAL
!     GPTA1  2.5.6
!     EST    2.3.56
!     SLT    2.3.51
 
!     ISTART GIVES 1ST POSITION OF HC OR REMFLUX VALUES
!     ROD IS IN ELEMENT COORDINATES, AS ARE TUBE,CONROD,BAR
 
!     X1 = 0.   X2 = X
!     AREA OF ROD NEEDED TO COMPUTE VOL
!     VOL  = LENGTH  * A
!     AREA OF TUBE CONPUTED WT OUTS.DIA.
 
!     OPEN FILE EST FOR ELEMENT DATA
 
!     INTEGRAL OVER VOL OF (GRADIENT SHAPE FUNC. TIMES GNU TIMES HC)
 
!     Z(1)  1ST POSITION OF LOAD
!     NELEMS = NO OF ELEMENTS
!     INCR   = MAX NO OF WORDS FOR AN ELEMENT OF THE ES T TABLE
!     NE(1 AND2) = ELEMENT NAME
 
 
 INTEGER, INTENT(IN)                      :: eltype
 INTEGER, INTENT(IN)                      :: istart
 INTEGER, INTENT(IN)                      :: itype
 INTEGER, INTENT(IN)                      :: ncount
 INTEGER, INTENT(IN)                      :: ido
 INTEGER, INTENT(IN)                      :: iwords
 INTEGER, INTENT(IN)                      :: nbdys
 INTEGER, INTENT(IN OUT)                  :: all
 INTEGER, INTENT(IN)                      :: nelout
 LOGICAL :: onlyc
 INTEGER :: estwds,outpt,sysbuf, scr6
 DIMENSION       xn(2),xload(2),nsil(2),iz(1),nam(2),  &
     necpt(200),NAME(2),hcx(2),hcy(2),hcz(2),  &
     zi(3),dndx(2),dndy(2),dndz(2),buf(50),ibuf(50),  &
     xlacc(3),xi(2),w(2),sc(5),isc(5)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ksystm(64)
 COMMON /zzzzzz/ z(1)
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /emecpt/ ecpt(200)
 COMMON /matin / matid,inflag,eltemp
 COMMON /hmtout/ xmat
 EQUIVALENCE     (ibuf(1)   ,buf(1)), (sc(1)  ,isc(1)  ),  &
     (ksystm(1 ),sysbuf), (ksystm(2),outpt ),  &
     (ksystm(56),ithrml), (ecpt(1),necpt(1)),  &
     (z(1)      ,iz(1) ), (nsil(1),necpt(2))
 DATA     nam  / 4H  EM,4H1D  /
 DATA     scr6 / 306          /
 
!     FROM EST GET ALL NECESSARY ELEMENT INFO
 
 onlyc  = .false.
 ng     = 3
 isc(1) = necpt(1)
 isc(2) = 1
 xi(1)  = -.5773502692
 xi(2)  = -xi(1)
 w(1)   = 1.
 w(2)   = 1.
 idx    = (eltype-1)*incr
 estwds = NE(idx+12)
 ngrids = NE(idx+10)
 NAME(1)= NE(idx+ 1)
 NAME(2)= NE(idx+ 2)
 
!     CHECK TO SEE IF THIS ELEMENT CONTAINS A GRID POINT ON A PERMBDY
!     CARD. IF SO, OR IF NO PERMBDY CARD EXISTS, COMPUTE LOADS FOR THE
!     ELEMENT IF NOT, COMPUTE HC CENTROIDAL VALUE ONLY. (ONLYC=.TRUE.)
!     THE PERMBDY SILS START AT Z(ISTART-NBDYS-1)
 
 IF (nbdys == 0) GO TO 20
 
 DO  i = 1,ngrids
   DO  j = 1,nbdys
     IF (nsil(i) == iz(istart-nbdys-nelout+j-1)) GO TO 20
   END DO
 END DO
 
!     ELEMENT HAS NO GRIDS ON PERMBDY
 
 onlyc = .true.
 ng    = 1
 20 CONTINUE
 IF (onlyc .AND. itype == 24) RETURN
 
!     IF ONLYC=TRUE, CHECK TO SEE IF THE ELEMENT HAD AN ELFORCE REQUEST.
!     IF SO, CONTINUE. IF NOT, JUST WRITE ZEROS TO HCCEN,SCR6) AND
!     RETURN.
 
 IF (.NOT.onlyc) GO TO 40
 IF (all == 1) GO TO 40
 IF (nelout == 0) GO TO 70
 
 DO  i = 1,nelout
   IF (necpt(1) == iz(istart-nelout+i-1)) GO TO 40
 END DO
 GO TO 70
 40 CONTINUE
 
!     1ST CHECK FOR ZERO LOAD
 
 IF (itype /= 20 .AND. itype /= 24) GO TO 80
 h1 = 0.
 h2 = 0.
 h3 = 0.
 DO  i = 1,2
   isub = istart + 3*nsil(i) - 3
   IF (itype == 24) isub = istart + 3*ncount - 3
   h1 = h1 + ABS(z(isub  ))
   h2 = h2 + ABS(z(isub+1))
   h3 = h3 + ABS(z(isub+2))
   IF (itype == 24) EXIT
 END DO
 60 hl = h1 + h2 + h3
 IF (hl    /= 0.) GO TO 80
 IF (itype == 24) RETURN
 
!     ALL ZEROS-WRITE ON SCR6
 
 70 sc(3) = 0.
 sc(4) = 0.
 sc(5) = 0.
 CALL WRITE (scr6,sc,5,0)
 RETURN
 
 80 CONTINUE
 
!     ROD ELTYPE 1
!     TUBE       3
!     CONROD    10
!     BAR       34
!     OTHERWISE  GET OUT
!     ONED  SOLVES LOAD DUE TO MAGNETIC FILED
 
 pi    = 3.14159
 inflag= 1
 IF (eltype /= 1 .AND. eltype /= 10) GO TO 90
 mid   = 4
 itemp = 17
 ix1   = 10
 ix2   = 14
 iy1   = 11
 iy2   = 15
 iz1   = 12
 iz2   = 16
 iar   = 5
 GO TO 110
 90 IF (eltype /= 3) GO TO 100
 mid   = 4
 itemp = 16
 ix1   = 9
 ix2   = 13
 iy1   = 10
 iy2   = 14
 iz1   = 11
 iz2   = 15
 
!     COMPUTE AREA
 
 dia   = ecpt(5)
 th    = ecpt(6)
 rad   = dia - 2.*th
 arrod = pi*((dia/2)**2 - (rad/2.)**2)
 GO TO 110
 100 IF (eltype /= 34) GO TO 300
 mid   = 16
 itemp = 42
 ix1   = 35
 ix2   = 39
 iy1   = 36
 iy2   = 40
 iz1   = 37
 iz2   = 41
 iar   = 17
 110 IF (onlyc) GO TO 120
 xl    = ecpt(ix2) - ecpt(ix1)
 yl    = ecpt(iy2) - ecpt(iy1)
 zl    = ecpt(iz2) - ecpt(iz1)
 xlen  = SQRT(xl**2 + yl**2 + zl**2)
 xn(1) = -1./xlen
 xn(2) =  1./xlen
 IF (eltype /= 3) arrod = ecpt(iar)
 eltemp= ecpt (itemp)
 matid = necpt(mid)
 
!     ARROD = AREA OF CROSS SECTION OF ROD
 
 IF (itype /= 24) CALL hmat (necpt(1))
 gnu = xmat
 IF (itype == 24) gnu = 1.
 
!     HC   FROM Z(ISTART)
!     YIELDS X COORD OF HC FOR GRID PT DEFINED BY NSIL
 
 vol = arrod*xlen
 
!     COMPUTE BASIC TO LOCAL TRANSFORMATION
 
 zi(1) = xl/xlen
 zi(2) = yl/xlen
 zi(3) = zl/xlen
 
!     PARTIALS OF N W.R.T X-GLOBAL,Y-GLOBAL,Z-GLOBAL
 
 dndx(1) = -zi(1)/xlen
 dndy(1) = -zi(2)/xlen
 dndz(1) = -zi(3)/xlen
 dndx(2) = -dndx(1)
 dndy(2) = -dndy(1)
 dndz(2) = -dndz(1)
 const   = .5*gnu*vol
 IF (itype == 24)GO TO 250
 120 CONTINUE
 jtype   = itype - 19
 xlacc(1)= 0.
 xlacc(2)= 0.
 xlacc(3)= 0.
 
!     LOOP OVER INTEGRATION POINTS-ASSUME CUBIC VARIATION. SO NEED 2
!     INTEGRATION POINTS + CENTROID
 
 DO  npts = 1,ng
   IF (npts /= ng) GO TO 130
   xx = .5*(ecpt(ix1) + ecpt(ix2))
   yy = .5*(ecpt(iy1) + ecpt(iy2))
   zz = .5*(ecpt(iz1) + ecpt(iz2))
   
!     AVERAGE SPCFLD
   
   xlx  = .5
   xlxp = .5
   GO TO 140
   
!     COMPUTE LOCAL COORDINATE OF SAMPLING POINT
   
   130 xlocal = .5*xlen*(1.+xi(npts))
   xlx  = xlocal/xlen
   xlxp = 1. - xlx
   
!     COMPUTE BASIC COORDS FOR XLOCAL
   
   xx   = xlxp*ecpt(ix1) + xlx*ecpt(ix2)
   yy   = xlxp*ecpt(iy1) + xlx*ecpt(iy2)
   zz   = xlxp*ecpt(iz1) + xlx*ecpt(iz2)
   140 ahcx = 0.
   ahcy = 0.
   ahcz = 0.
   
!     COMPUTE HC AT THIS POINT DUE TO ALL LOADS OF THIS TYPE
   
   DO  ijk = 1,ido
     IF (itype == 20) GO TO 160
     isub = istart + (ijk-1)*iwords - 1
     DO  i = 1,iwords
       buf(i) = z(isub+i)
     END DO
     SELECT CASE ( jtype )
       CASE (    1)
         GO TO 160
       CASE (    2)
         GO TO 180
       CASE (    3)
         GO TO 190
       CASE (    4)
         GO TO 200
     END SELECT
     
!     SPCFLD
     
     160 DO  i = 1,2
       is = istart + 3*nsil(i) -3
       hcx(i) = z(is  )
       hcy(i) = z(is+1)
       hcz(i) = z(is+2)
     END DO
     
!     INTERPOLATE GRID VALUES TO INTEGRATION POINT
     
     hc1 = xlxp*hcx(1) + xlx*hcx(2)
     hc2 = xlxp*hcy(1) + xlx*hcy(2)
     hc3 = xlxp*hcz(1) + xlx*hcz(2)
     GO TO 210
     180 CALL axloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 210
     190 CALL geloop (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     GO TO 210
     200 CALL dipole (buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     210 ahcx  = ahcx + hc1
     ahcy  = ahcy + hc2
     ahcz  = ahcz + hc3
   END DO
   IF (npts /= ng) GO TO 230
   sc(3) = ahcx
   sc(4) = ahcy
   sc(5) = ahcz
   CALL WRITE (scr6,sc,5,0)
   IF (onlyc) RETURN
   CYCLE
   
!     WE HAVE HC AT THIS INTEGRATION POINT. MULT. BY WEIGHT AND
!     ACCUMULATE
   
   230 xlacc(1) = xlacc(1) + ahcx*w(npts)
   xlacc(2) = xlacc(2) + ahcy*w(npts)
   xlacc(3) = xlacc(3) + ahcz*w(npts)
 END DO
 
!     MULT. BY CONST AND GRAD N TO GET LOADS
 
 xload(1) = const*(dndx(1)*xlacc(1) + dndy(1)*xlacc(2) + dndz(1)*xlacc(3))
 xload(2) = const*(dndx(2)*xlacc(1) + dndy(2)*xlacc(2) + dndz(2)*xlacc(3))
 GO TO 260
 
!     REMFLUX
 
 250 is   = istart + 3*ncount - 3
 ahcx = z(is  )
 ahcy = z(is+1)
 ahcz = z(is+2)
 
 xload(1) = gnu*vol*(dndx(1)*ahcx + dndy(1)*ahcy + dndz(1)*ahcz)
 xload(2) = gnu*vol*(dndx(2)*ahcx + dndy(2)*ahcy + dndz(2)*ahcz)
 260 DO  i = 1,2
   j = nsil(i)
   
!     IF PERMBDY EXISTS AND IF GRID IS NOT ON IT, IGNORE ITS LOAD
   
   IF (nbdys == 0) GO TO 280
   DO  k = 1,nbdys
     IF (j /= iz(istart-nbdys-nelout+k-1)) CYCLE
     GO TO 280
   END DO
   GO TO 290
   280 CONTINUE
   290 z(j) = z(j) - xload(i)
 END DO
 RETURN
 
 300 WRITE  (outpt,310) ufm
 310 FORMAT (a23,', ELEMENT TYPE ',2A4,' WAS USED IN AN E AND M ', 'PROBLEM.')
 CALL mesage (-37,0,nam)
 RETURN
END SUBROUTINE em1d
