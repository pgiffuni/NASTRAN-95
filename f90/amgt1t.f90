SUBROUTINE amgt1t (nline,nl,acpt,nstns,c3t,c4t)
     
    !     GENERATE CONSTANTS C3T AND C4T FOR
    !     STREAMLINE NL OF SWEPT TURBOPROP BLADE.
 
 
    INTEGER, INTENT(IN OUT)                  :: nline
    INTEGER, INTENT(IN)                      :: nl
    INTEGER, INTENT(IN OUT)                  :: acpt
    INTEGER, INTENT(IN OUT)                  :: nstns
    REAL, INTENT(OUT)                        :: c3t
    REAL, INTENT(OUT)                        :: c4t
    REAL :: l1,l2,l3
    INTEGER :: FILE,NAME(2)
    DIMENSION  pn(3),p1(3),fn(3),f1(3),s1(3),sn(3),DATA(3)
    DATA FILE/ 102 /, NAME / 4HAMGT,4H1T  /
 
    !     INPUT VARIABLES -
    !     NLINE      TOTAL NO. OF STREAMLINES
    !     NL         PRESENT STEAMLINE
    !     ACPT       SCRATCH UNIT WITH BASIC COORDINATES OF NODES
    !     NSTNS      TOTAL NO. OF STATIONS
 
    !     OUTPUT VARIABLES -
    !     C3T        CONSTANTS USED BY SUB. AND SUP. AERODYNAMIC ROUTINES
    !     C4T        USED IN DEFINING DATA BLOCK AJJ
 
    !     LOCAL VARIABLES -
    !     PN         COORDINATES TRAILING EDGE PREVIOUS STREAMLINE
    !     P1         COORDINATES LEADING EDGE PREVIOUS STREAMLINE
    !     FN         COORDINATES TRAILING EDGE NEXT STREAMLINE
    !     F1         COORDINATES LEADING EDGE NEXT STREAMLINE
    !     S1         COORDINATES OF LEADING EDGE OF CURRENT STREAMLINE
    !     SN         COORDINATES OF TRAILING EDGE OF CURRENT STREAMLINE
 
    !     EXTRACT LEADING COORDINATES OF CURRENT STREAMLINE
 
    CALL fread (acpt,DATA,3,0)
    DO  i = 1,3
        s1(i) = DATA(i)
    END DO
 
    !     SKIP TO TRAILING EDGE COORDINATES OF CURRENT STREAMLINE
 
    nskip = (2-nstns)*3
    CALL READ (*905,*900,acpt,DATA,nskip,0,mm)
    CALL fread (acpt,DATA,3,0)
    DO  i = 1,3
        sn(i) = DATA(i)
    END DO
 
    !     EXTRACT COORDINATES FOR PREVIOUS--P-FOR FIRST STREAMLINE
 
    IF (nl /= 1) GO TO 10
    DO  i = 1,3
        pn(i) = sn(i)
        p1(i) = s1(i)
    END DO
 
    !     NOW COORDINATES FOR NEXT -F- FOR LAST STREAMLINE
 
10  IF (nl /= nline) GO TO 15
    DO  i = 1,3
        fn(i) = sn(i)
        f1(i) = s1(i)
    END DO
    GO TO 50
 
    !     NOW COORDINATES FOR NEXT -F- FOR ALL OTHER STREAMLINES
 
    !     SKIP FIRST 10 WORDS OF NEXT STREAMLINE
 
15  CALL READ (*905,*900,acpt,DATA,-10,0,mm)
    CALL fread (acpt,DATA,3,0)
    f1(1) = DATA(1)
    f1(2) = DATA(2)
    f1(3) = DATA(3)
 
    !     COMPUTE SKIP TO TRAILING EDGE COORDINATES
 
    nskip = (2-nstns)*3
    CALL READ (*905,*900,acpt,DATA,nskip,0,mm)
    CALL fread (acpt,DATA,3,0)
    fn(1) = DATA(1)
    fn(2) = DATA(2)
    fn(3) = DATA(3)
50  a1 = sn(1) - s1(1)
    b1 = sn(2) - s1(2)
    c1 = sn(3) - s1(3)
 
    a2 = fn(1) - p1(1)
    b2 = fn(2) - p1(2)
    c2 = fn(3) - p1(3)
 
    a3 = pn(1) - f1(1)
    b3 = pn(2) - f1(2)
    c3 = pn(3) - f1(3)
 
    a4 = b2*c1 - b1*c2
    b4 = c2*a1 - c1*a2
    c4 = a2*b1 - a1*b2
 
    a5 = b1*c3 - b3*c1
    b5 = c1*a3 - c3*a1
    c5 = a1*b3 - a3*b1
 
    l1 = SQRT(a1**2 + b1**2 + c1**2)
    l2 = SQRT(a4**2 + b4**2 + c4**2)
    l3 = SQRT(a5**2 + b5**2 + c5**2)
 
    a6 = 0.5*(a4/l2 + a5/l3)
    b6 = 0.5*(b4/l2 + b5/l3)
    c6 = 0.5*(c4/l2 + c5/l3)
 
    a7 = (b1*c6 - b6*c1)/l1
    b7 = (c1*a6 - c6*a1)/l1
    c7 = (a1*b6 - a6*b1)/l1
 
    a8 = f1(1) - p1(1)
    b8 = f1(2) - p1(2)
    c8 = f1(3) - p1(3)
 
    a9 = fn(1) - pn(1)
    b9 = fn(2) - pn(2)
    c9 = fn(3) - pn(3)
 
    w1 = a7*a8 + b7*b8 + c7*c8
    w2 = a7*a9 + b7*b9 + c7*c9
 
    c3t = (w2-w1)/(2.0*l1)
    c4t = w1/2.0
 
    IF (nl == nline) RETURN
 
    !     RETURN TO START OF RECORD
 
    CALL bckrec (acpt)
 
    !     COMPUTE SKIP TO NEXT STREAMLINE AT EXIT FROM THIS ROUTINE
 
    nskip = -6 - (10+3*nstns)*nl
    CALL READ (*905,*900,acpt,DATA,nskip,0,mm)
 
    !     SET PREVIOUS COORDINATES -P- TO PRESENT STREAMLINE COORDINATES
 
    DO  i = 1,3
        pn(i) = sn(i)
        p1(i) = s1(i)
    END DO
    RETURN
 
    !     E-O-R    ENCOUNTERED
900 ip1 = -3
    GO TO 999
 
    !     E-O-F    ENCOUNTERED
 
905 ip1 = -2
999 CALL mesage (ip1,FILE,NAME)

    RETURN
END SUBROUTINE amgt1t
