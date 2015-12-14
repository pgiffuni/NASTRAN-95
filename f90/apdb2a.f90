SUBROUTINE apdb2a (nline,nl,scr1,nstns,m1,s1,sn,tblt,tblr)
     
!     GENERATE BASIC TO LOCAL TRANSFORMATION MATRIX FOR
!     STREAMLINE NL OF SWEPT TURBOPROP BLADE.
 
 
 INTEGER, INTENT(IN OUT)                  :: nline
 INTEGER, INTENT(IN)                      :: nl
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(IN)                      :: nstns
 REAL, INTENT(IN)                         :: m1
 REAL, INTENT(IN)                         :: s1(3)
 REAL, INTENT(IN)                         :: sn(3)
 REAL, INTENT(OUT)                        :: tblt(3)
 REAL, INTENT(OUT)                        :: tblr(3)
 REAL :: l1,l2,l3
 
 INTEGER :: FILE,NAME(2)
 
 DIMENSION pn(3),p1(3),fn(3),f1(3)
 DIMENSION DATA(7)
 
 DATA FILE/301/,NAME  /4HAPDB,4H2A  /
 
!---------------------------------------------------------------------
!     INPUT VARIABLES--
!     NLINE      TOTAL NO. OF STREAMLINES
!     NL         PRESENT STEAMLINE
!     SCR1       SCRATCH UNIT WITH BASIC COORDINATES OF NODES
!     NSTNS      TOTAL NO. OF STATIONS
!     M1         SIGN BASED ON ROTATION OF BLADE
!     S1         COORDINATES OF LEADING EDGE OF CURRENT STREAMLINE
!     SN         COORDINATES OF TRAILING EDGE OF CURRENT STREAMLINE
 
!     OUTPUT VARIABLES--
!     TBLT      BASIC TO LOCAL TRANSFORMATION FOR TRANSLATION
!     TBLR      BASIC TO LOCAL TRANSFORMATION FOR ROTATION
 
!     LOCAL VARIABLES--
!     PN        COORDINATES TRAILING EDGE PREVIOUS STREAMLINE
!     P1        COORDINATES LEADING EDGE PREVIOUS STREAMLINE
!     FN        COORDINATES TRAILING EDGE NEXT STREAMLINE
!     F1        COORDINATES LEADING EDGE NEXT STREAMLINE
!---------------------------------------------------------------------
!     EXTRACT COORDINATES FOR PREVIOUS--P-FOR FIRST STREAMLINE
!---------------------------------------------------------------------
 IF(nl /= 1)GO TO 10
 DO  i=1,3
   pn(i)=sn(i)
   p1(i)=s1(i)
 END DO
!----------------------------------------------------------------------
!     NOW COORDINATES FOR NEXT--F-FOR LAST STREAMLINE
!----------------------------------------------------------------------
 10   IF(nl /= nline)GO TO 15
 DO  i=1,3
   fn(i)=sn(i)
   f1(i)=s1(i)
 END DO
 GO TO 50
!----------------------------------------------------------------------
!     NOW COORDINATES FOR NEXT--F-FOR ALL OTHER STREAMLINES
!---------------------------------------------------------------------
 15   CALL fread(scr1,DATA,7,0)
 f1(1)=DATA(5)
 f1(2)=DATA(6)
 f1(3)=DATA(7)
!----------------------------------------------------------------------
!    COMPUTE SKIP TO TRAILING EDGE COORDINATES
!-----------------------------------------------------------------------
 nskip=(2-nstns)*7
 CALL READ(*905,*900,scr1,DATA,nskip,0,mm)
 CALL fread(scr1,DATA,7,0)
 fn(1)=DATA(5)
 fn(2)=DATA(6)
 fn(3)=DATA(7)
!----------------------------------------------------------------------
!     RETURN TO START OF RECORD
!----------------------------------------------------------------------
 CALL bckrec(scr1)
!---------------------------------------------------------------------
!     COMPUTE SKIP TO ORIGINAL LOCATION AT ENTRY TO THIS ROUTINE
!---------------------------------------------------------------------
 nskip=-nstns*nl*7
 CALL READ(*905,*900,scr1,DATA,nskip,0,mm)
 50   a1=sn(1)-s1(1)
 b1=sn(2)-s1(2)
 c1=sn(3)-s1(3)
 
 a2=fn(1)-p1(1)
 b2=fn(2)-p1(2)
 c2=fn(3)-p1(3)
 
 a3=pn(1)-f1(1)
 b3=pn(2)-f1(2)
 c3=pn(3)-f1(3)
 
 a4=b2*c1-b1*c2
 b4=c2*a1-c1*a2
 c4=a2*b1-a1*b2
 
 a5=b1*c3-b3*c1
 b5=c1*a3-c3*a1
 c5=a1*b3-a3*b1
 
 l1=SQRT(a1**2+b1**2+c1**2)
 l2=SQRT(a4**2+b4**2+c4**2)
 l3=SQRT(a5**2+b5**2+c5**2)
 
 a6=0.5 *(a4/l2  +  a5/l3)
 b6=0.5 *(b4/l2  +  b5/l3)
 c6=0.5 *(c4/l2  +  c5/l3)
!---------------------------------------------------------------------
!     BASIC TO LOCAL TRANSFORMATION FOR TRANSLATION
!---------------------------------------------------------------------
 tblt(1)= a6*m1
 tblt(2)= b6*m1
 tblt(3)= c6*m1
!----------------------------------------------------------------------
!     BASIC TO LOCAL TRANSFORMATION FOR ROTATION
!---------------------------------------------------------------------
 tblr(1)= -m1*a1/l1
 tblr(2)= -m1*b1/l1
 tblr(3)= -m1*c1/l1
 IF(nl == nline)RETURN
!---------------------------------------------------------------------
!     SET PREVIOUS COORDINATES--P- TO PRESENT STREAMLINE
!---------------------------------------------------------------------
 DO  i=1,3
   pn(i)=sn(i)
   p1(i)=s1(i)
 END DO
 RETURN
!     E-O-R    ENCOUNTERED
 900  ip1 = -3
 GO TO 999
!     E-O-F    ENCOUNTERED
 905  ip1 = -2
 999  CALL mesage(ip1,FILE,NAME)
 RETURN
END SUBROUTINE apdb2a
