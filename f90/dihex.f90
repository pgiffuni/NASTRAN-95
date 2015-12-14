SUBROUTINE dihex (TYPE)
     
!     THIS ROUTINE PROCESSES IHEX1, IHEX2, AND IHEX3 ELEMENT DATA TO
!     PRODUCE THE DIFFERENTIAL STIFFNESS MATRIX
 
!           TYPE = 1    IHEX1
!           TYPE = 2    IHEX2
!           TYPE = 3    IHEX3
 
!     THE EST ENTRIES ARE
 
!     NAME  ----------INDEX----------   DESCRIPTION
!            IHEX1    IHEX2    IHEX3
 
!     EID        1        1        1    ELEMENT ID NO.
!     SIL      2-9     2-21     2-33    SCALAR INDEX LIST
!     MID       10       22       34    MATERIAL ID NO.
!     CID       11       23       35    MATERIAL COORD. SYSTEM ID NO.
!     NIP       12       24       36    NO. INTEGRATION POINTS PER EDGE
!     MAXAR     13       25       37    MAX ASPECT RATIO
!     ALFA      14       26       38    MAX ANGLE FOR NORMALS
!     BETA      15       27       39    MAX ANGLE FOR MIDSIDE POINTS
!     BGPDT  16-47   28-107   40-167    BASIC GRID POINT DATA
!     GPT    48-55  108-127  168-199    GRID POINT TEMPERATURES
!     DEF       56      128      200    NOT USED
!     GPTLD  57-64  129-148  201-232    GRID POINT TEMPERATURE LOADS
!     UGV    65-88  149-208  233-328    GLOBAL DISPLACEMENT VECTOR
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL :: anis       ,rect       ,tdep       ,diag       , mtdep
 INTEGER :: heat       ,eid        ,sil(1)     ,ib(46)     ,  &
      jz(1)      ,cid        ,iest(1)    ,  &
     bcord      ,bgpdt      ,gpt        ,nc(8)      ,  &
     ufm(6)     ,elno(3)    ,iwork(1)   ,otpt
 REAL :: nu         ,maxar      ,dmaxar(3)  ,dalfa(3)   ,  &
     dbeta(2)   ,evec(3,12) ,work(66)   ,vn(3,2)
 DOUBLE PRECISION :: sk         ,sv         ,z(1)       ,jacob(3,3) ,  &
     detj       ,s(4)       ,h(4)       ,gauss(8)   ,  &
     sfact      ,part(3,3)  ,e1         ,e2         ,  &
     e3         ,tf(3,3)    ,tk(3,3)    ,sig(6)     ,  &
     sx         ,sy         ,sz         ,sxy        ,  &
     syz        ,szx        ,str(18)    ,c(3,3)
 DOUBLE PRECISION :: gmat(36)    ,store(18)  ,dalpha(6)
 DIMENSION        igrid(128) ,grid(128)
 COMMON /matin/   mid        ,inflag     ,temp
 COMMON /matout/  e          ,g          ,nu         ,rho        ,  &
     talpha     ,tref       ,cdamp      ,SPACE(18)  , mtdep
 COMMON /matiso/  bufm6(46)
 COMMON /ds1aaa/  npvt       ,dum6(6)    ,i6x6k      ,dum12(12)  ,  &
     jmax       ,dum2(2)    ,nrowsc
 COMMON /zzzzzz/  zs(1)
 COMMON /ds1adp/  isil(32)   ,sk(6,6)    ,work       ,str
 COMMON /ds1aet/  est(328)
 COMMON /system/  sysbuf,otpt,nogo,sys(6),mtemp
 EQUIVALENCE      (z(1),jz(1),zs(1))    ,(eid,est(1),iest(1))   ,  &
     (sil(1),est(2))       ,(work(1),iwork(1))     ,  &
     (sig(1),sx)           ,(sig(2),sy)            ,  &
     (sig(3),sz)           ,(sig(4),sxy)           ,  &
     (sig(5),syz)          ,(sig(6),szx)
 EQUIVALENCE      (work(1),evec(1,1))   ,(work(37),vn(1,1))     ,  &
     (work(43),nc(1))
 EQUIVALENCE      (work(1),jacob(1,1))  ,(work(19),h(1))        ,  &
     (work(27),s(1))       ,(work(35),part(1,1))   ,  &
     (work(53),sig(1))     ,(work(1),c(1,1))
 EQUIVALENCE      (work(1),tf(1,1))     ,(work(35),tk(1,1))
 EQUIVALENCE      (ib(1),bufm6(1))      ,(grid(1),igrid(1))
 DATA    kgg   /  101 /, mgg   / -1     /
 DATA    dmaxar,  dalfa, dbeta / 5.0    ,10.0      ,15.0        ,  &
     45.0    ,45.0      ,45.0        , 45.0      ,45.0        /
 DATA    dtor  ,  gauss  /    0.017453292519943E0, 0.577350269189626D0,  &
     0.555555555555556D0, 0.774596669241483D0,  &
     0.888888888888889D0, 0.347854845137454D0,  &
     0.861136311594053D0, 0.652145154862546D0,  &
     0.339981043584856D0/
 DATA    ufm          /4H0***,4H use,4HR fa,4HTAL ,4HMESS,4HAGE /
 DATA    ihex , elno  /4HIHEX,4H ele,4HMENT,4H no./
 DATA    nerr1        / 2141 /
 
 heat = 0
 
!     FOR DOUBLE PRECISION, OPEN CORE POINTERS MUST BE MODIFIED
 
 izs = 1 + 2*(i6x6k + jmax*nrowsc)
 nzs = izs + 10655
 iz  = izs/2 + 1
 nz  = nzs/2 + 1
 iprec = 2
 
!     ALLOCATE LARGE ARRAYS IN OPEN CORE
 
 ngp = 12*TYPE - 4
 DO  i = 1,ngp
   IF (sil(i) == npvt) GO TO 7
 END DO
 nogo= 1
 7 igp = i
 IF (heat == 1) GO TO 30
 ngg = 3*ngp
 IF (kgg <= 0) GO TO 10
 ik  = iz + 3*ngg
 nk  = ik - 1 + (ngg+1)*ngg/2
 GO TO 20
 10 ik  = iz
 nk  = ik + 3*ngg - 1
 im  = nk + 1
 nm  =(ngp+1)*ngp/2 + nk
 GO TO 40
 20 nm  = nk
 IF (mgg <= 0) GO TO 40
 im  = nk + 1
 nm  = nk + (ngp+1)*ngp/2
 GO TO 40
 30 ik  = iz + 17
 nk  = ik - 1 + ngp**2
 im  = nk + 1
 nm  = im - 1 + ngp**2
 ngg = ngp
 40 in  = nm + 1
 ig  = in + ngp
 ix  = ig + 3*ngp
 nd  = nm + 9*ngp
 id  = nd + 1
 nd  = id + ngg - 1
 IF (nd <= nz) GO TO 100
 WRITE (otpt,7100) ufm,nerr1,ihex,TYPE,elno,eid
 nogo = 1
 
!     OPEN CORE MAP
!     =============
 
!     DOUBLE PRECISION  Z(1)
!     COMMON  /ZZZZZZ/  Z
 
!     NGG = ORDER OF ELEMENT MATRIX
 
!     INDEX      STIFFNESS             MASS                HEAT
!                AND MASS              ONLY              TRANSFER
 
!     IZ    NGG BY 3 PARTITION  NGG BY 3 PARTITION  FOUR WORD COORDINATE
!           OF MATRIX           OF MATRIX           VECTOR.  INPUT TO
!                                                   TRANSD
 
!     IZ+2                                          TRANSFORMED THERMAL
!                                                   CONDUCTANCE MATRIX
 
!     IT                                            MATERIAL TRANSFOR-
!                                                   MATION MATRIX
 
!     IK    SYMMETRIC HALF OF   SAME AS IZ          FULL CONDUCTANCE
!           STIFFNESS
 
!     IM    SYMMETRIC HALF OF   SYMMETRIC HALF OF   FULL CAPACITANCE
!           MASS                MASS
 
!     IN    --------------------SHAPE FUNCTIONS-------------------------
 
!     IG    --------------------D(SHAPE)/D(GREEK)-----------------------
 
!     IX    --------------------D(SHAPE)/D(BASIC XYZ)-------------------
 
!     ID    DISPLACEMENT
!           VECTOR IN BASIC
!           COORDINATES
 
 
!     CHECK GEOMETRY.  THE FOLLOWING CHECKS ARE MADE
!           1.  ASPECT RATIO
!           2.  ANGLES BETWEEN NORMALS OF SUB-TRIANGLES ON EACH FACE
!           3.  ANGLES BETWEEN VECTORS BETWEEN POINTS ALONG EACH EDGE
!           4.  REVERSE SEQUENCING
!           5.  DUPLICATE COORDINATE VALUES
 
!     FETCH EPT DATA, COMPUTE EST POINTERS
 
 100 mid  = 10 + 12*(TYPE-1)
 cid  = iest(mid+1)
 nip  = iest(mid+2)
 maxar= est(mid+3)
 alfa = est(mid+4)
 beta = est(mid+5)
 bgpdt= mid + 6
 gpt  = bgpdt + 4*ngp
 mid  = iest(mid)
 IF (nip < 2 .OR. nip > 4) nip = TYPE/2 + 2
 IF (maxar <= 0.0) maxar = dmaxar(TYPE)
 IF (alfa  < 0.0) alfa  = dalfa(TYPE)
 IF (beta  < 0.0) beta  = dbeta(TYPE-1)
 alfa = COS(dtor*alfa)
 beta = COS(dtor*beta)
 
!     TRANSFORM DISPLACEMENT VECTOR TO BASIC COORDINATES
 
 DO  i = 1,ngp
   m  = bgpdt + 4*i - 4
   n  = id + 3*i - 3
   j  = gpt + 2*ngp + 3*(i-1) + 1
   IF (iest(m) == 0) GO TO 102
   CALL transd (est(m) ,z(iz))
   DO  l = 1,3
     z(ik+l-1) = DBLE(est(j+l-1)*0.25)
   END DO
   CALL gmmatd (z(iz),3,3,0,z(ik),3,1,0,z(n))
   CYCLE
   102 DO  l = 1,3
     z(n+l-1) = DBLE(est(j+l-1)*0.25)
   END DO
 END DO
 
!     REARRANGE BGPDT
 
 DO  i = 1,ngp
   igrid(i) = iest(bgpdt+4*i-4)
 END DO
 bcord = gpt - 3
 DO  i = 2,ngp
   DO  j = 1,3
     k = bgpdt + 4*(ngp-i) + 4 - j
     bcord = bcord - 1
     est(bcord) = est(k)
   END DO
 END DO
 DO  i = 2,ngp
   iest(bgpdt+i-1) = igrid(i)
 END DO
 
!     INITIALIZE FOR NUMERICAL INTEGRATION
 
 
!     ABSCISSAE AND WEIGHT COEFFICIENTS FOR GAUSSIAN QUADRATURE
 
 i = nip - 1
 SELECT CASE ( i )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
 END SELECT
 510 h(1) = 1.0
 s(1) = gauss(1)
 h(2) = 1.0
 s(2) =-gauss(1)
 GO TO 540
 520 h(1) = gauss(2)
 s(1) = gauss(3)
 h(2) = gauss(4)
 s(2) = 0.0
 h(3) = gauss(2)
 s(3) =-gauss(3)
 GO TO 540
 530 h(1) = gauss(5)
 s(1) = gauss(6)
 h(2) = gauss(7)
 s(2) = gauss(8)
 h(3) = gauss(7)
 s(3) =-gauss(8)
 h(4) = gauss(5)
 s(4) =-gauss(6)
 
!     GENERATE TABLE OF EQUIVALENTS IN SIL ARRAY SO MATRIX WILL BE
!     ORDERED ACCORDING TO INCREASING SIL NUMBERS
 
 540 DO  i = 1,ngp
   isil(i) = sil(i)
   sil(i) = i
 END DO
 
!     NOW SIL(I) = PARTITION NUMBER OF ELEMENT GRID POINT I
 
!     ZERO OUT OPEN CORE FOR MATRIX SUMMATION
 
 DO  i = ik,nm
   z(i) = 0.0
 END DO
 
!     FETCH MATERIAL PROPERTIES
 
!     THIS SECTION OF CODE MUST BE UPDATED WHEN GENERAL ANISOTROPIC
!     MATERIAL IS ADDED.
 
!     TEST FOR ANISOTROPIC MATERIAL
 
 anis = .false.
 
!     TEST FOR RECTANGULAR COORDINATE SYSTEM IN WHICH THE ANISOTROPIC
!     MATERIAL IS DEFINED
 
 rect = .true.
 
!     CHECK FOR TEMPERATURE DEPENDENCE
 
 tdep = .true.
 DO  i = 2,ngp
   IF (est(gpt) /= est(gpt+i-1)) GO TO 630
 END DO
 tdep = .false.
 630 temp = est(gpt)
 inflag = 10
 CALL mat (eid)
 IF (.NOT.mtdep) tdep = .false.
 IF (ib(46) == 6) anis = .true.
 tref = bufm6(44)
 
 IF (kgg <= 0) GO TO 1000
 
!     IF ISOTROPIC, TEMPERATURE INDEPENDENT MATERIAL, COMPUTE CONSTANTS
 
 IF (tdep) GO TO 1000
 IF (anis) GO TO 800
 IF (ib(46) /= 0) GO TO 640
 WRITE (otpt,7300) ufm,mid,eid
 nogo = 1
 RETURN
 
 640 e1 = bufm6(1)
 e2 = bufm6(2)
 e3 = bufm6(22)
 talpha = bufm6(38)
 GO TO 1000
 
!     IF MATERIAL IS ANISOTROPIC, DEFINED IN A RECTANGULAR COORDINATE
!     SYSTEM, AND NOT TEMPERATURE DEPENDENT, TRANSFORM IT TO BASIC
!     SYSTEM
 
 800 DO  ijk = 1,36
   gmat(ijk) = bufm6(ijk)
 END DO
 
!     CODE TO TRANSFORM GENERAL ANISOTROPIC MATERIAL PROPERTIES TO
!     BASIC COORDINATE SYSTEM MUST BE ADDED HERE.
 
!     ALL SET TO BEGIN INTEGRATION LOOPS.  DO IT.
 
 1000 DO  i = 1,nip
   DO  j = 1,nip
     DO  k = 1,nip
       
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE
       
       CALL ihexsd (TYPE,z(in),z(ig),jacob,detj,eid,s(i),s(j),s(k),  &
           est(bcord))
       IF (detj /= 0.0) GO TO 1010
       
!     BAD ELEMENT IF FALL HERE.  JACOBIAN MATRIX WAS SINGULAR.
       
       nogo = 1
       RETURN
       
       1010 sfact = h(i)*h(j)*h(k)*detj
       IF (kgg <= 0) GO TO 1015
       
!     STIFFNESS
       
!     COMPUTE STRAIN-DISPLACEMENT RELATIONS
       
!     MUST REVERSE CALLING ORDER SINCE MATRICES ARE STORED BY COLUMNS
       
       CALL gmmatd (z(ig),ngp,3,0,jacob,3,3,0,z(ix))
       
!     IF MATERIAL IS TEMPERATURE DEPENDENT, MUST COMPUTE TEMPERATURE
!     AT THIS INTEGRATION POINT AND FETCH MATERIAL PROPERTIES AGAIN
       
       1015 IF (.NOT.tdep) GO TO 1030
       temp = 0.0
       DO  l = 1,ngp
         temp = temp + z(in+l-1)*est(gpt+l-1)
       END DO
       CALL mat (eid)
       IF (kgg <= 0) GO TO 1100
       IF (anis) GO TO 1040
       IF (ib(46) /= 0 ) GO TO 1025
       WRITE (otpt,7300) ufm,mid,eid
       nogo = 1
       RETURN
       
       1025 e1 = bufm6(1)
       e2 = bufm6(2)
       e3 = bufm6(22)
       talpha = bufm6(38)
       GO TO 1100
       1030 IF (kgg <= 0) GO TO 1100
       
!     IF MATERIAL IS ANISOTROPIC AND NOT DEFINED IN RECTANGULAR COOR-
!     DINATE SYSTEM, MUST TRANSFORM TO BASIC COORDINATE SYSTEM AT THIS
!     INTEGRATION POINT
       
!     THIS CODE MUST BE COMPLETED WHEN GENERAL ANISOTROPIC MATERIAL IS
!     ADDED
       
       IF (.NOT.anis) GO TO 1100
       IF (rect) GO TO 1100
       1040 CONTINUE
       
!     INSERT GLOBAL TO BASIC TRANSFORMATION OPERATIONS HERE FOR
!     ANISOTROPIC MATERIAL MATRIX
       
       DO  ijk = 1,36
         gmat(ijk) = bufm6(ijk)
       END DO
       
!     MATERIAL HAS BEEN EVALUATED FOR THIS INTEGRATION POINT WHEN
!     FALL HERE.
       
       1100 CONTINUE
       
!     COMPUTE STRESSES FOR DIFFERENTIAL STIFFNESS MATRIX
       
!     THERMAL EFFECTS
       
       IF (iest(gpt+ngp+1) == -1) GO TO 1120
       
!     COMPUTE LOADING TEMPERATURE AT THIS POINT
       
       temp = 0.0
       DO  l = 1,ngp
         temp = temp + z(in+l-1)*est(gpt+ngp+l)
       END DO
       temp = 0.25*(temp-tref)
       IF (anis) GO TO 1115
       sig(1) =-talpha*(e1+2.0*e2)*temp
       sig(2) = sig(1)
       sig(3) = sig(1)
       sig(4) = 0.0
       sig(5) = 0.0
       sig(6) = 0.0
       GO TO 1140
       
!     ANISOTROPIC
       
       1115 DO  ijk = 1,6
         dalpha(ijk) = bufm6(ijk+37)
       END DO
       
       CALL gmmatd (gmat(1),6,6,0,dalpha(1),6,1,0,sig)
       DO  ijk = 1,6
         sig(ijk) = -sig(ijk)*temp
       END DO
       GO TO 1140
       1120 DO  l = 1,6
         sig(l) = 0.0
       END DO
       
!     DISPLACEMENT EFFECTS, COMPUTE STRESS MATRIX AND MULTIPLY BY DISPL.
       
       1140 str(12) = 0.0
       str(13) = 0.0
       str(17) = 0.0
       DO  l = 1,ngp
         ii = ix + 3*l - 4
         IF (anis) GO TO 1141
         str( 1) = e1*z(ii+1)
         str( 2) = e2*z(ii+2)
         str( 3) = e2*z(ii+3)
         str( 4) = e2*z(ii+1)
         str( 5) = e1*z(ii+2)
         str( 6) = e2*z(ii+3)
         str( 7) = e2*z(ii+1)
         str( 8) = e2*z(ii+2)
         str( 9) = e1*z(ii+3)
         str(10) = e3*z(ii+2)
         str(11) = e3*z(ii+1)
         str(14) = e3*z(ii+3)
         str(15) = e3*z(ii+2)
         str(16) = e3*z(ii+3)
         str(18) = e3*z(ii+1)
         GO TO 1145
         
!     ANISOTROPIC
         
         1141 DO  ijk = 1,18
           store(ijk) = 0.d0
         END DO
         store( 1) = z(ii+1)
         store( 5) = z(ii+2)
         store( 9) = z(ii+3)
         store(10) = z(ii+2)
         store(11) = z(ii+1)
         store(14) = z(ii+3)
         store(15) = z(ii+2)
         store(16) = z(ii+3)
         store(18) = z(ii+1)
         
         CALL gmmatd (gmat(1),6,6,0,store(1),6,3,0,str(1))
         
         1145 CONTINUE
         CALL gmmatd (str,6,3,-2,z(id+3*l-3),3,1,0,sig)
       END DO
       sv = sx
       sx = sx + sy
       sy = sy + sz
       sz = sz + sv
       
!     NOW BEGIN LOOPS OVER GRID POINTS ALONG ROWS AND COLUMNS
       
       DO  n = 1,ngp
         DO  m = n,ngp
           IF (n == igp .OR. m == igp) GO TO 1170
           CYCLE
           1170 CONTINUE
           
!     COMPUTE PARTITION FOR POINTWISE ROW M AND COLUMN N
           
           IF (kgg <= 0) GO TO 1300
           IF (.NOT.anis ) GO TO 1200
           
!     MUST ADD CODE TO COMPUTE THE CONTRIBUTION TO THE STIFFNESS MATRIX
!     FOR ANISOTROPIC MATERIAL HERE
!     =================================================================
           
           1200 IF (sil(m) >= sil(n)) GO TO 1210
           
!     MUST COMPUTE TRANSPOSE OF THIS PARTITION FOR SUMMATION IN ELEMENT
!     MATRIX
           
           mz = ix + (n-1)*3
           nz = ix + (m-1)*3
           GO TO 1220
           1210 mz = ix + (m-1)*3
           nz = ix + (n-1)*3
           
!     DIFFERENTIAL STIFFNESS
           
           1220 DO  l = 1,3
             DO  inc = 1,3
               c(l,inc)  = z(mz+inc-1)*z(nz+l-1)
             END DO
           END DO
           part(1,1) = sx*c(2,2) + syz*(c(2,3) + c(3,2)) + sz*c(3,3)
           part(2,2) = sy*c(3,3) + szx*(c(3,1) + c(1,3)) + sx*c(1,1)
           part(3,3) = sz*c(1,1) + sxy*(c(1,2) + c(2,1)) + sy*c(2,2)
           part(2,1) =-sx*c(2,1) + sxy*c(3,3) - syz*c(1,3) - szx*c(2,3)
           part(3,1) =-sz*c(3,1) - sxy*c(3,2) - syz*c(2,1) + szx*c(2,2)
           part(1,2) =-sx*c(1,2) + sxy*c(3,3) - syz*c(3,1) - szx*c(3,2)
           part(3,2) =-sy*c(3,2) - sxy*c(3,1) + syz*c(1,1) - szx*c(1,2)
           part(1,3) =-sz*c(1,3) - sxy*c(2,3) - syz*c(1,2) + szx*c(2,2)
           part(2,3) =-sy*c(2,3) - sxy*c(1,3) + syz*c(1,1) - szx*c(2,1)
           
!     ADD STIFFNESS PARTITION TO ELEMENT MATRIX
           
!     COMPUTE INDEX INTO OPEN CORE WHERE PART(1,1) IS TO BE ADDED.
           
           IF (sil(m)-sil(n) < 0.0) THEN
             GO TO  1230
           ELSE IF (sil(m)-sil(n) == 0.0) THEN
             GO TO  1240
           ELSE
             GO TO  1250
           END IF
           1230 mz = sil(n)
           nz = sil(m)
           diag = .false.
           GO TO 1260
           1240 mz = sil(m)
           nz = sil(n)
           diag = .true.
           GO TO 1260
           1250 mz = sil(m)
           nz = sil(n)
           diag = .false.
           
!     COLUMN NUMBER
           
           1260 l = (nz-1)*3 + 1
           
!     INCREMENT BETWEEN COLUMNS
           
           inc = ngg - l
           
!     FIRST WORD OF COLUMN
           
           l = ik + ((l-1)*l)/2 + (inc+1)*(l-1)
           
!     WORD IN COLUMN FOR THIS ROW
           
           l = l + 3*(mz-nz)
           
!     ADD PARTITION
           
           DO  nz = 1,3
             DO  mz = 1,3
               IF (diag .AND. mz < nz) CYCLE
               z(l+mz-1) = z(l+mz-1) + part(mz,nz)*sfact
             END DO
             l   = l + inc
             inc = inc - 1
           END DO
           1300 IF (mgg <= 0) CYCLE
         END DO
       END DO
     END DO
   END DO
 END DO
 
!     END OF INTEGRATION LOOPS
 
 
!     LOOK FOR NON-BASIC COORDINATE SYSTEM
 
 DO  i = 1,ngp
   IF (iest(bgpdt+i-1) /= 0) GO TO 2005
 END DO
 GO TO 2061
 
!     RESTORE GRID POINT DATA TO ORIGINAL FORM FOR DOING TRANSFORM
!     TO GLOBAL COORDINATES
 
!     FIRST, TRANSFER IT TO OPEN CORE AT IN
 
 2005 k = (in-1)*2 + 1
 j = ngp*4
 DO  i = 1,j
   grid(i) = est(bgpdt+i-1)
 END DO
 
!     NOW MOVE IT BACK AND REARRANGE IT
 
 DO  i = 1,ngp
   iest(bgpdt+4*i-4) = igrid(i)
   DO  j = 1,3
     est(bgpdt+4*i-4+j) = grid(ngp+3*i+j-3)
   END DO
 END DO
 
!     FETCH GLOBAL TO BASIC TRANSFORMATION MATRICES
 
 DO  i = 1,ngp
   j = in + (i-1)*9
   CALL transd (est(bgpdt+4*i-4),z(j))
 END DO
 
!     TRANSFORM STIFFNESS TO GLOBAL COORDINATES
 
 DO  i = 1,ngp
   
!     COLUMN INDICES
   
   k   = (i-1)*3 + 1
   inc = ngg - k + 1
   l   = ik + ((k-1)*k)/2 + inc*(k-1)
   m   = l + inc
   n   = m + inc - 1
   
!     TRANSFORMATION MATRIX INDEX
   
   nz  = in + (i-1)*9
   
!     TERMS ON DIAGONAL PARTITION
   
   CALL tktztk (tk,z,nz,l,m,n)
   
!     OFF-DIAGONAL PARTITIONS
   
   l   = l + 3
   m   = m + 2
   n   = n + 1
   irp = i + 1
   IF (irp > ngp) CYCLE
   DO  j = irp,ngp
     nz  = in + 9*(j-1)
     DO  k = 1,3
       tk(k,1) = 0.0
       tk(k,2) = 0.0
       tk(k,3) = 0.0
       DO  mz = 1,3
         tk(k,1) = tk(k,1) + z(l+mz-1)*z(nz+3*mz+k-4)
         tk(k,2) = tk(k,2) + z(m+mz-1)*z(nz+3*mz+k-4)
         tk(k,3) = tk(k,3) + z(n+mz-1)*z(nz+3*mz+k-4)
       END DO
     END DO
     mz = in + 9*(i-1)
     DO  k = 1,3
       z(l+k-1) = 0.0
       z(m+k-1) = 0.0
       z(n+k-1) = 0.0
       DO  ii = 1,3
         z(l+k-1) = z(l+k-1) + tk(k,ii)*z(mz+3*ii-3)
         z(m+k-1) = z(m+k-1) + tk(k,ii)*z(mz+3*ii-2)
         z(n+k-1) = z(n+k-1) + tk(k,ii)*z(mz+3*ii-1)
       END DO
     END DO
     l = l + 3
     m = m + 3
     n = n + 3
   END DO
 END DO
 
!     BUILD STIFFNESS PARTITIONS AND PASS TO EMGOUT
 
 2061 DO  i = 1,36
   sk(i,1) = 0.0D0
 END DO
 i = igp
 DO  j = 1,3
   
!     COLUMN NUMBER
   
   k = (i-1)*3 + j
   
!     NUMBER OF TERMS TO FETCH TO COMPLETE THIS COLUMN IN PARTITION
   
   l = k - 1
   IF (l == 0) GO TO 2075
   
!     FETCH TERMS AND LOAD INTO J-TH COLUMN OF PARTITION
   
   n   = ik + l
   inc = ngg - 1
   DO  m = 1,l
     z(iz+ngg*j-ngg+m-1) = z(n)
     n   = n + inc
     inc = inc - 1
   END DO
   
!     FILL OUT PARTITION WITH COLUMNS OF STIFFNESS MATRIX
   
!     COMPUTE INDEX IN OPEN CORE OF FIRST TERM OF COLUMN K
   
   2075 n = ik + ((k-1)*k)/2 + (ngg-k+1)*(k-1)
   
!     INSERT THIS COLUMN IN PARTITION
   
   DO  m = k,ngg
     z(iz+ngg*j-ngg+m-1) = z(n)
     n = n + 1
   END DO
 END DO
 DO  i = 1,ngp
   j = (i-1)*3 + iz - 1
   DO  m = 1,3
     DO  n = 1,3
       k = j + n + (m-1)*ngg
       sk(n,m) = z(k)
     END DO
   END DO
   CALL ds1b (sk,isil(i))
 END DO
 
!     ALL DONE, NO ERRORS
 
 RETURN
 
 
 7100 FORMAT (6A4,i4,2H, ,a4,i1,3A4,i9,' INSUFFICIENT CORE TO COMPUTE',  &
     ' ELEMENT MATRIX')
 7300 FORMAT (6A4,'4005. AN ILLEGAL VALUE OF -NU- HAS BEEN SPECIFIED ',  &
     'UNDER MATERIAL ID =',i10,17H for element id =,i10)
 
END SUBROUTINE dihex
