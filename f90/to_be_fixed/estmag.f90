SUBROUTINE estmag (hest,estfld,mpt,dit,geom1,iany,kcount)
     
!     CREATE SCRATCH FILE ESTFLD WHICH WILL BE USED TO COMPUTE TOTAL
!     MAGNETIC FIELD. READ EST AND CREATE SIMILAR RECORDS CONTAINING
!     ELTYPE,EID,NUMBER OF SILS,SILS,3 X 3 MATERIALS MATRIX, AND 3 X 3
!     TRANSFORMATION MATRIX TO BRING HM BACK TO BASIC COORD. SYSTEM
!     FROM ELEMENT SYSTEM,OUTPUT COORD. SYSTEM ID, AND BASIC COORDS.
!     OF STRESS POINT(USUALLY AVERAGE OF GRID COORDS.)
 
 
 INTEGER, INTENT(IN)                      :: hest
 INTEGER, INTENT(IN)                      :: estfld
 INTEGER, INTENT(IN)                      :: mpt
 INTEGER, INTENT(IN)                      :: dit
 INTEGER, INTENT(IN)                      :: geom1
 INTEGER, INTENT(OUT)                     :: iany
 INTEGER, INTENT(OUT)                     :: kcount
 INTEGER :: sysbuf,eltype,pointr(6,20),mcb(7),  &
     buf1,buf2,iz(1),nam(2),estwds,frstgd,otpe, bfield(2),oldeid,buf3,FILE
 DIMENSION       dn(8),xm(32),coord(3),kount(2),NAME(2),v12(3),  &
     v13(3),xi(3),xj(3),xk(3),ecpt(200),iecpt(200), e(9),g(9)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /gpta1 / nelems,last,incr,NE(1)
 COMMON /system/ sysbuf,otpe
 COMMON /zzzzzz/ z(1)
 COMMON /hmatdd/ iihmat,nnhmat,mptfil,iditfl
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /hmtout/ xmat(6)
 EQUIVALENCE     (z(1),iz(1)),(ecpt(1),iecpt(1))
 DATA    nam   / 4HESTM,4HAG  /
 DATA    bfield/ 3101,31      /
 DATA    dn    / 4*-.25, 4*.5 /
 
!     ITYPE ITH MID ISYS1 ITEMP DIM
 
 DATA pointr/   1,   0,  4,   9,   17,   1, 3,   0,  4,   8,   16,   1,  &
     6,   5,  6,  15,   27,   2, 9,   5,  6,   9,   21,   2,  &
     10,   0,  4,   9,   17,   1, 16,   6,  7,  10,   26,   2,  &
     17,   5,  6,   9,   21,   2, 18,   6,  7,  10,   26,   2,  &
     19,   6,  7,  16,   32,   2, 34,   0, 16,  34,   42,   1,  &
     36,   5,  6,   7,   19,   2, 37,   6,  7,   8,   24,   2,  &
     39,   0,  2,   7,   23,   3, 40,   0,  2,   9,   33,   3,  &
     41,   0,  2,  11,   43,   3, 42,   0,  2,  11,   43,   3,  &
     65,   0, 10,  16,   48,   3, 66,   0, 22,  28,  108,   3,  &
     67,   0, 34,  40,  168,   3, 80,  11, 12,  14,   46,   2     /
 
 lcore = korsz(z)
 buf1  = lcore - sysbuf + 1
 buf2  = buf1  - sysbuf
 buf3  = buf2  - sysbuf - 1
 lcore = buf3  - 1
 IF (lcore <= 0) GO TO 1008
 
!     COUNT NUMBER OF TERMS IN A COLUMN OF HCCEN
 
 kcount = 0
 
!     SET UP MATERIALS
 
 iihmat = 0
 nnhmat = lcore
 mptfil = mpt
 iditfl = dit
 CALL prehma (z)
 nextz  = nnhmat + 1
 
 CALL gopen (hest,z(buf1),0)
 CALL gopen (estfld,z(buf2),1)
 
!     READ IN ANY BFIELD CARDS
 
 iany   = 0
 iall   = 1
 ifield = 0
 idefid = 0
 nfield = 0
 FILE   = geom1
 CALL preloc (*1001,z(buf3),geom1)
 CALL locate (*3,z(buf3),bfield,idex)
 iany   = 1
 CALL READ (*1002,*1,geom1,z(nextz+1),lcore-nextz,0,iwords)
 GO TO 1008
 1 nfield = iwords/2
 IF (nfield /= 1 .OR. iz(nextz+2) /= -1) GO TO 2
 ifield = iz(nextz+1)
 idefid = ifield
 GO TO 3
 
!     BFIELD ARE NOT THE SAME FOR EVERY ELEMENT
 
 2 iall = 0
 3 CALL CLOSE (geom1,1)
 
!     CHECK FOR ALL SO THAT CSTM WONT BE OPENED
 
 IF (nfield == 0) GO TO 8
 DO  i = 1,iwords,2
   IF (iz(nextz+i) /= 0) GO TO 5
 END DO
 iany = 0
 iall = 1
 ifield = 0
 5 CONTINUE
 
!     CHECK FOR A DEFAULT ID
 
 DO  i = 2,iwords,2
   IF (iz(nextz+i) == -1) GO TO 7
 END DO
 GO TO 8
 7 idefid = iz(nextz+i-1)
 8 CONTINUE
 FILE = hest
 
 10 CALL READ (*120,*1003,hest,eltype,1,0,iflag)
 CALL WRITE (estfld,eltype,1,0)
 oldeid = 0
 icount = 0
 idx    = (eltype-1)*incr
 estwds = NE(idx+12)
 ngrids = NE(idx+10)
 frstgd = 2
 IF (eltype >= 39 .AND. eltype <= 42) frstgd = 3
 NAME(1) = NE(idx+1)
 NAME(2) = NE(idx+2)
 
!     PICK UP MATERIAL ID, START OF BGPDT DATA, AND DIMENSIONALITY OF
!     ELEMENT
 
 DO  i = 1,20
   jel = i
   IF (eltype-pointr(1,i) < 0.0) THEN
     GO TO   500
   ELSE IF (eltype-pointr(1,i) == 0.0) THEN
     GO TO    30
   ELSE
     GO TO    20
   END IF
 END DO
 GO TO 500
 
 30 ith   = pointr(2,jel)
 mid   = pointr(3,jel)
 isys1 = pointr(4,jel)
 isys2 = isys1 + 4
 isys3 = isys2 + 4
 
!     FOR IS2D8, USE 4TH POINT FOR GEOMETRY SINCE THAT IS WHAT WE USE
!     FOR IS2D8 ELSEWHERE
 
 IF (eltype == 80) isys3 = isys3 + 4
 itemp = pointr(5,jel)
 idim  = pointr(6,jel)
 
 40 CALL READ (*1002,*110,hest,ecpt,estwds,0,iflag)
 IF (eltype < 65) kcount = kcount + 3
 IF (eltype == 65) kcount = kcount + 27
 IF (eltype == 66 .OR. eltype == 67) kcount = kcount + 63
 IF (eltype == 80) kcount = kcount + 27
 
!     FIND BFIELD FOR THIS ELEMENT
 
 IF (iall == 1) GO TO 47
 DO  i = 2,iwords,2
   IF (iecpt(1) == iz(nextz+i)) GO TO 46
 END DO
 ifield = idefid
 GO TO 47
 46 ifield = iz(nextz+i-1)
 47 CONTINUE
 
!     WRITE EID, SILS
 
 CALL WRITE (estfld,iecpt(1),1,0)
 CALL WRITE (estfld,ngrids,1,0)
 CALL WRITE (estfld,iecpt(frstgd),ngrids,0)
 
!     FETCH MATERIALS
 
 matid = iecpt(mid)
 sinth = 0.
 costh = 0.
!***
!    ASSUME HERE THAT FOR ISOPARAMETRICS WE HAVE TEMPERATURE-INDEPENDENT
!    MATERIALS IN THIS MAGNETICS PROBLEM
!***
 eltemp = ecpt(itemp)
 inflag = 3
 CALL hmat (iecpt(1))
 g(1) = xmat(1)
 g(2) = xmat(2)
 g(3) = xmat(3)
 g(4) = xmat(2)
 g(5) = xmat(4)
 g(6) = xmat(5)
 g(7) = xmat(3)
 g(8) = xmat(5)
 g(9) = xmat(6)
 
!     NOW CREATE TRANSFORMATION MATRIX FROM LOACL COORDS TO BASIC
!     DETERMINE DIMENSIONALITY OF ELEMENT
 
 SELECT CASE ( idim )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 70
   CASE (    3)
     GO TO 90
 END SELECT
 
!     ONE-DIMENSIONAL-DETERMINE THE LOCAL X-AXIS(IN BASIC COORDS)
 
 50 DO  i = 1,3
   v12(i) = ecpt(isys2+i) - ecpt(isys1+i)
 END DO
 xlen   = SQRT(v12(1)**2 + v12(2)**2 + v12(3)**2)
 DO  i = 1,3
   v12(i) = v12(i)/xlen
 END DO
 
 DO  i = 1,9
   e(i) = 0.
 END DO
 e(1) = v12(1)
 e(4) = v12(2)
 e(7) = v12(3)
 GO TO 100
 
!     TWO-DIMENSIONAL  WE WILL USE ONLY GRIDS 1,2,3 ASSUMING A PLANAR
!     OR NEARLY PLANAR ELEMENT FOR QUADS
 
 70 DO  i = 1,3
   v12(i) = ecpt(isys2+i) - ecpt(isys1+i)
   v13(i) = ecpt(isys3+i) - ecpt(isys1+i)
 END DO
 xlen  = SQRT(v12(1)**2 + v12(2)**2 + v12(3)**2)
 DO  i = 1,3
   xi(i) = v12(i)/xlen
 END DO
 
 xk(1) = xi(2)*v13(3) - xi(3)*v13(2)
 xk(2) = xi(3)*v13(1) - xi(1)*v13(3)
 xk(3) = xi(1)*v13(2) - xi(2)*v13(1)
 xlen  = SQRT(xk(1)**2 + xk(2)**2 + xk(3)**2)
 DO  i = 1,3
   xk(i) = xk(i)/xlen
 END DO
 
 xj(1) = xk(2)*xi(3) - xk(3)*xi(2)
 xj(2) = xk(3)*xi(1) - xk(1)*xi(3)
 xj(3) = xk(1)*xi(2) - xk(2)*xi(1)
 xlen  = SQRT(xj(1)**2 + xj(2)**2 + xj(3)**2)
 DO  i = 1,3
   xj(i) = xj(i)/xlen
 END DO
 
 DO  i = 1,3
   e(3*i-2) = xi(i)
   e(3*i-1) = xj(i)
   e(3*i  ) = xk(i)
 END DO
 
!     CHECK ON MATERIALS AS IN EMRING
 
 angle = ecpt(ith)*0.017453293
 IF (xmat(3) == 0. .AND. xmat(5) == 0.) GO TO 87
 GO TO 100
 87 IF (ABS(angle) <= .0001) GO TO 100
 DO  i = 1,9
   g(i) = 0.
 END DO
 s    = SIN(angle)
 c    = COS(angle)
 csq  = c*c
 ssq  = s*s
 cs   = c*s
 x2   = 2.*cs*xmat(2)
 g(1) = csq*xmat(1) - x2 + ssq*xmat(4)
 g(2) = cs*(xmat(1) - xmat(4)) + (csq-ssq)*xmat(2)
 g(5) = ssq*xmat(1) + x2 + csq*xmat(4)
 g(4) = g(2)
 g(9) = xmat(6)
 
 IF (eltype /= 36 .AND. eltype /= 37) GO TO 100
 
!     SINCE MAT5 INFO FOR TRAPRG,TRIARG MUST BE GIVEN IN X-Y TERMS,
!     RE-ORDER THE 3 X 3 , INTERCHANGING Y AND Z
 
 temp = g(5)
 g(5) = g(9)
 g(9) = temp
 temp = g(2)
 g(2) = g(3)
 g(3) = temp
 g(4) = g(2)
 g(7) = g(3)
 GO TO 100
 
!     THREE-DIMENSIONAL-NO ELEMENT COORDINATE SYSTEM-EVERYTHING IS
!     OUTPUT IN BASIC- SO E IS IDENTITY
 
 90 DO  i = 1,9
   e(i) = 0.
 END DO
 e(1) = 1.
 e(5) = 1.
 e(9) = 1.
 
 100 CALL WRITE (estfld,g,9,0)
 CALL WRITE (estfld,e,9,0)
 IF (eltype >= 65 .AND. eltype <= 67) GO TO 104
 
!     COMPUTE THE AVERAGE COORDINATES OF THE GRID POINTS OF THE ELEMENT
!     FOR USE IN NON-RECTANGULAR COORDIANTE SYSTEMS. THIS POINT IS NOT
!     NECESSARILY THE CENTROID,BUT ANY POINT WILL DO FOR CONSTANT STRAIN
!     ELEMENTS AND THIS IS CONVENIENT
 
 IF (eltype /= 80) GO TO 1013
 
!     FOR IS2D8 USE SHAPE FUNCTION
 
 DO  i = 1,3
   coord(i) = 0.
   DO  j = 1,8
     jsub = isys1 + 4*(j-1)
     coord(i) = coord(i) + dn(j)*ecpt(jsub+i)
   END DO
 END DO
 GO TO 108
 1013 CONTINUE
 DO  i = 1,3
   coord(i) = 0.
   DO  j = 1,ngrids
     jsub = isys1 + 4*(j-1)
     coord(i) = coord(i) + ecpt(jsub+i)
   END DO
 END DO
 DO  i = 1,3
   coord(i) = coord(i)/FLOAT(ngrids)
 END DO
 GO TO 108
 
!     ISOPARAMETRICS-PICK UP COORDS. OF APPLICABLE POINT. FOR THE LAST
!     POINT, GO TO THE PREVIOUS METHOD
 
 104 IF (iecpt(1) == oldeid) GO TO 105
 oldeid = iecpt(1)
 105 icount = icount + 1
 IF (eltype == 65 .AND. icount < 9) GO TO 106
 IF (eltype > 65 .AND. icount < 21) GO TO 106
 
!     CENTROIDAL POINT-COMPUTE COORDS BASED ON XI=ETA=ZETA=0
 
 icount = 0
 oldeid = 0
 IF (eltype /= 65) GO TO 1051
 DO  i = 1,8
   xm(i) = .125
 END DO
 GO TO 1057
 1051 IF (eltype /= 66) GO TO 1054
 DO  i = 1,20
   xm(i) = .25
 END DO
 DO  i = 1,7,2
   xm(i) =-.25
   xm(i+12) =-.25
 END DO
 GO TO 1057
 1054 con1 =  9./64.
 con2 =-19./64.
 DO  i = 1,32
   xm(i) = con1
 END DO
 DO  i = 1,10,3
   xm(i) = con2
   xm(i+20) = con2
 END DO
 
 1057 DO  i = 1,3
   coord(i) = 0.
   DO  j = 1,ngrids
     jsub = isys1 + 4*(j-1)
     coord(i) = coord(i) + ecpt(jsub+i)*xm(j)
   END DO
 END DO
 GO TO 109
 106 IF (eltype == 67 .AND. icount < 21) GO TO 1071
 jsub = isys1 + 4*(icount-1)
 DO  i = 1,3
   coord(i) = ecpt(jsub+i)
 END DO
 IF (icount > 1) GO TO 109
 GO TO 108
 
!     FOR IHEX3, MUST GET PROPER COORDINATES IF NOT THE LAST POINT
 
 1071 IF (icount >= 9 .AND. icount <= 12) GO TO 1072
 IF ((icount/2)*2 == icount) GO TO 1073
 
!     CORNERS
 
 IF (icount == 1 .OR. icount == 13) jcount =-1
 iadd = 0
 IF (icount >= 13) iadd = 8
 jcount = jcount + 1
 num = 1
 kount(1) = icount + jcount + iadd
 GO TO 1075
 
!     MIDSIDES
 
 1072 kadd = 4
 jco  = 3
 GO TO 1074
 1073 kadd = 1
 IF (icount == 2 .OR. icount == 14) jco = -1
 1074 iadd = 0
 IF (icount >= 14) iadd = 8
 jco  = jco + 1
 num  = 2
 kount(1) = icount + jco + iadd
 kount(2) = kount(1) + kadd
 1075 DO  i = 1,3
   coord(i) = 0.
   DO  j = 1,num
     jsub = isys1 + 4*(kount(j)-1)
     coord(i) = coord(i) + ecpt(jsub+i)
   END DO
 END DO
 DO  i = 1,3
   coord(i) = coord(i)/FLOAT(num)
 END DO
 IF (icount > 1) GO TO 109
 
!     WRITE OUT CID AND COORDINATES
 
 108 CALL WRITE (estfld,ifield,1,0)
 109 CALL WRITE (estfld,coord,3,0)
 
!     FOR ISOPARAMETRICS, GET COORDS OF NEXT POINT, OTHERWISE,
!     GO BACK FOR ANOTHER ELEMENT OF THIS TYPE
 
 IF (oldeid == 0) GO TO 40
 GO TO 105
 
 
!     GET ANOTHER ELEMENT TYPE
 
 110 CALL WRITE (estfld,0,0,1)
 GO TO 10
 
!     DONE
 
 120 CALL CLOSE (estfld,1)
 CALL CLOSE (hest,1)
 mcb(1) = hest
 CALL rdtrl (mcb)
 mcb(1) = estfld
 CALL wrttrl (mcb)
 RETURN
 
!     FATAL ERRORS
 
 500 WRITE  (otpe,501) ufm,NAME
 501 FORMAT (a23,', ELEMENT TYPE ',2A4,' NOT ALLOWED IN ESTMAG')
 CALL mesage (-61,0,0)
 
 1001 n = -1
 GO TO 1010
 1002 n = -2
 GO TO 1010
 1003 n = -3
 GO TO 1010
 1008 n = -8
 FILE = 0
 1010 CALL mesage (n,FILE,nam)
 RETURN
END SUBROUTINE estmag
