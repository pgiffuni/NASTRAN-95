SUBROUTINE ihex(temps,pg,TYPE)
     
!     ELEMENT THERMAL LOAD GENERATOR FOR ISOPARAMETRIC SOLID ELEMENTS
 
!     TYPE = 1     CIHEX1
!     TYPE = 2     CIHEX2
!     TYPE = 3     CIHEX3
 
!***********************************************************************
!           THE EST ENTRIES ARE
 
!     NAME  ---------INDEX---------   DESCRIPTION
!            IHEX1   IHEX2   IHEX3
 
!     EID        1       1       1    ELEMENT ID NO.
!     SIL      2-9    2-21    2-33    SCALAR INDEX LIST
!     MID       10      22      34    MATERIAL ID NO.
!     CID       11      23      35    MATERIAL COORD. SYSTEM ID NO.
!     NIP       12      24      36    NO. INTEGRATION POINTS PER EDGE
!     MAXAR     13      25      37    MAX ASPECT RATIO
!     ALFA      14      26      38    MAX ANGLE FOR NORMALS
!     BETA      15      27      39    MAX ANGLE FOR MIDSIDE POINTS
!     BGPDT  16-47  28-107  40-167    BASIC GRID POINT DATA
!     GPT    48-55 108-127 168-199    GRID POINT TEMPERATURES
!***********************************************************************
 
 
 REAL, INTENT(IN OUT)                     :: temps(1)
 REAL, INTENT(OUT)                        :: pg(1)
 INTEGER, INTENT(IN)                      :: TYPE
 LOGICAL :: tdep     ,mtdep      ,anis       ,rect
 
 INTEGER :: otpt       ,eid        ,iest(1)    ,bgpdt      ,  &
     bcord    ,gpt        ,jz(32)     ,sil        ,cid        , ufm(6)
 INTEGER :: ib(46)
 
 DOUBLE PRECISION :: shp        ,dshp       ,jacob      ,detj  &
     ,     s          ,sfact      ,a(6)      ,e1         ,e2  &
     ,     e3         ,parg(96)   ,cn(3,32)   ,temp       ,eltemp  &
     ,     alpvec     ,gmat(36)   ,gauss(8)   ,dalpha(6)
 
 REAL :: psgl(96)
 
 COMMON/trimex/ est(200)
 COMMON/matin/ mid      ,inflag     ,eltemp
 COMMON/matout/   se          ,g          ,snu        ,rho        ,  &
     talpha,tref,cdamp,SPACE(18), mtdep
 COMMON/matiso/ bufm6(46)
 COMMON/system/ sysbuf  ,otpt       ,sys1(7)          ,mtemp      ,  &
     sys2(45),heat
 
 COMMON/ssgwrk/    shp(32)    ,dshp(3,32) ,jacob(3,3) ,s(4) ,     h(4)
 
 EQUIVALENCE (eid,est(1),iest(1)),(jz(1),shp(1))
 EQUIVALENCE (psgl(1),parg(1))
 EQUIVALENCE (ib(1),bufm6(1))
 
 DATA  gauss/      0.577350269189626D0    ,0.555555555555556D0  &
     ,                 0.774596669241483D0    ,0.888888888888889D0  &
     ,                 0.347854845137454D0    ,0.861136311594053D0  &
     ,                 0.652145154862546D0    ,0.339981043584856D0/
 DATA ufm /4H0***,4H use,4HR fa,4HTAL ,4HMESS,4HAGE /
 
!*****
!     COMPUTE EST POINTERS
!*****
 ngp = 12*TYPE - 4
 mid = 10 + 12*(TYPE - 1)
 cid=iest(mid+1)
 nip=iest(mid+2)
 IF (nip < 2 .OR. nip > 4) nip=TYPE/2+2
 bgpdt = mid + 6
 gpt=bgpdt+4*ngp
 DO  i=1,ngp
   jz(i) = iest(bgpdt + 4*i - 4)
 END DO
 bcord=gpt-3
 DO   i=2,ngp
   DO   j=1,3
     k = bgpdt + 4*(ngp - i) + 4 - j
     bcord = bcord - 1
     est(bcord) = est(k)
   END DO
 END DO
 DO   i=2,ngp
   iest(bgpdt+i-1) = jz(i)
 END DO
 mid=iest(ngp+2)
 
!     ABSCISSAE AND WEIGHT COEFFICIENTS FOR GAUSSIAN QUADRATURE
 
 i=nip-1
 SELECT CASE ( i )
   CASE (    1)
     GO TO 131
   CASE (    2)
     GO TO 132
   CASE (    3)
     GO TO 133
 END SELECT
 131 h(1)=1.0
 s(1)=gauss(1)
 h(2)=h(1)
 s(2)=-s(1)
 GO TO 134
 132 h(1)=gauss(2)
 s(1)=gauss(3)
 h(2)=gauss(4)
 s(2)=0.0
 h(3)=h(1)
 s(3)=-s(1)
 GO TO 134
 133 h(1)=gauss(5)
 s(1)=gauss(6)
 h(2)=gauss(7)
 s(2)=gauss(8)
 h(3)=h(2)
 s(3)=-s(2)
 h(4)=h(1)
 s(4)=-s(1)
 134 CONTINUE
 
!=======================================================================
!     THIS SECTION OF CODE MUST BE UPDATED WHEN GENERAL ANISOTROPIC
!     MATERIAL IS ADDED
 
!     TEST FOR ANISOTROPIC MATERIAL
 
 anis = .false.
 inflag=10
 
!     TEST FOR RECTANGULAR COORDINATE SYSTEM IN WHICH THE ANISOTROPIC
!     MATERIAL IS DEFINED
 
 rect = .true.
!=======================================================================
 
!     FETCH MATERIAL AND SET TEMPERATURE DEPENDENCE FLAG
 
 tdep=.true.
 DO  i=2,ngp
   IF (est(gpt) /= est(gpt+i-1)) GO TO 150
 END DO
 tdep=.false.
 150 eltemp=est(gpt)
 CALL mat(eid)
 IF (.NOT. mtdep) tdep=.false.
 IF (ib(46) == 6) anis=.true.
 tref=bufm6(44)
!*****
!     IF ISOTROPIC TEMPERATURE INDEPENDENT MATERIAL, COMPUTE CONSTANTS
!*****
 IF (tdep) GO TO 800
 IF (anis) GO TO 700
 IF (ib(46) /= 0) GO TO 640
 CALL page2(2)
 WRITE(otpt,7300) ufm,mid,eid
 nogo = 1
 RETURN
 640 e1=bufm6(1)
 e2=bufm6(2)
 e3=bufm6(22)
 talpha=bufm6(38)
 GO TO 800
 
!=======================================================================
!     CODE TO TRANSFORM GENERAL ANISOTROPIC MATERIAL PROPERTIES TO
!     BASIC COORDINATE SYSTEM MUST BE ADDED HERE
!=======================================================================
 
 700 DO  ijk=1,36
   gmat(ijk)=bufm6(ijk)
 END DO
 800 ntlp = 3*ngp
 DO  i=1,ntlp
   parg(i) = 0.0
 END DO
!*****
!     BEGIN INTEGRATION LOOP NOW
!*****
 DO  i=1,nip
   DO  j=1,nip
     DO  k=1,nip
!*****
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE
!*****
       CALL ihexsd(TYPE,shp,dshp,jacob,detj,eid,s(i),s(j),s(k), est(bcord))
       IF (detj /= 0.0D0) GO TO 1010
       
!     JACOBIAN MATRIX WAS SINGULAR
       
       CALL mesage(-61,0,0)
!*****
!     COMPUTE PARTIAL DERIVATIVE OF SHAPE FUNCTIONS WITH RESPECT
!     TO BASIC COORDINATES
!*****
       1010 CALL gmmatd(dshp,ngp,3,0,jacob,3,3,0,cn)
!*****
!     COMPUTE LOADING TEMPERATURE AT THIS INTEGRATION POINT
!*****
       temp=0.0D0
       DO  l=1,ngp
         temp=temp+shp(l)*DBLE(temps(l))
       END DO
       temp=temp-DBLE(tref)
!*****
!     IF MATERIAL IS TEMPERATURE DEPENDENT, COMPUTE TEMPERATURE AT THIS
!     INTEGRATION POINT AND FETCH MATERIAL PROPERTIES
!*****
       IF(.NOT.tdep)  GO TO 1030
       eltemp=0.0D0
       DO   l=1,ngp
         eltemp=eltemp+shp(l)*DBLE(est(gpt+l-1))
       END DO
       CALL mat(eid)
       IF (anis) GO TO 1040
       IF (ib(46) /= 0) GO TO 1025
       CALL page2(2)
       WRITE(otpt,7300)  ufm,mid,eid
       nogo = 1
       RETURN
       1025 e1=bufm6(1)
       e2=bufm6(2)
       e3=bufm6(22)
       talpha=bufm6(38)
       GO TO 1100
!*****
!     IF MATERIAL IS ANISOTROPIC AND NOT DEFINED IN RECTANGULAR COOR-
!     DINATE SYSTEM, MUST TRANSFORM TO BASIC COORDINATE SYSTEM AT THIS
!     INTEGRATION POINT
!*****
       1030 IF(.NOT. anis)  GO TO 1100
       IF (rect) GO TO 1500
       1040 CONTINUE
       DO  ijk=1,36
         gmat(ijk)=bufm6(ijk)
       END DO
       
!=======================================================================
!     INSERT GLOBAL TO BASIC TRANSFORMATION OPERATIONS HERE FOR
!     ANISOTROPIC MATERIAL MATRIX
       GO TO 1500
!=======================================================================
!*****
!     COMPUTE CONTRIBUTION TO THERMAL LOAD VECTOR FOR ISOTROPIC MATERIAL
!*****
       1100 alpvec=DBLE(talpha)*(e1+2.0*e2)
       sfact=h(i)*h(j)*h(k)*detj*alpvec*temp
       l = 0
       DO   ii=1,ngp
         DO  jj=1,3
           l = l + 1
           parg(l) = sfact*cn(jj,ii)  + parg(l)
         END DO
       END DO
       CYCLE
!=======================================================================
       1500 CONTINUE
!     ADD LOAD COMPUTATIONS FOR ANISOTROPIC MATERIAL HERE
!=======================================================================
       
       sfact=h(i)*h(j)*h(k)*detj*temp
       DO  ijk=1,6
         dalpha(ijk)=bufm6(ijk+37)
       END DO
       
       CALL gmmatd(gmat,6,6,0,dalpha,6,1,0,a(1))
       l=0
       DO  ii=1,ngp
         l=l+1
         parg(l)=parg(l)+sfact*(cn(1,ii)*a(1)+cn(2,ii)*a(4)+cn(3,ii)*a(6))
         l=l+1
         parg(l)=parg(l)+sfact*(cn(2,ii)*a(2)+cn(1,ii)*a(4)+cn(3,ii)*a(5))
         l=l+1
         parg(l)=parg(l)+sfact*(cn(3,ii)*a(3)+cn(2,ii)*a(5)+cn(1,ii)*a(6))
       END DO
     END DO
   END DO
 END DO
 DO  i=1,ntlp
   psgl(i)=parg(i)
 END DO
!*****
!     INSERT THERMAL LOAD INTO GLOBAL LOAD VECTOR  (PG ARRAY)
!*****
 
 DO  i=1,ngp
   sil = iest(i+1)
   ibgp = bgpdt + i - 1
   IF (iest(ibgp) == 0) GO TO 2500
   CALL basglb(psgl(3*i-2),psgl(3*i-2),est(bcord+3*i-3),iest(ibgp))
   2500 DO  j=1,3
     pg(sil+j-1)=pg(sil+j-1)+psgl(3*i-3+j)
   END DO
 END DO
 
 
 7300 FORMAT(6A4,69H4005. an illegal value of -nu- has been specified un  &
     der material id =,i10,17H for element id =,i10)
 RETURN
END SUBROUTINE ihex
