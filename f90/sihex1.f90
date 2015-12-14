SUBROUTINE sihex1 (TYPE,strspt,nip)
     
!     PHASE 1 STRESS ROUTINE FOR IHEX1, IHEX2, AND IHEX3 ELEMENTS
 
!     TYPE = 1    IHEX1
!     TYPE = 2    IHEX2
!     TYPE = 3    IHEX3
 
!     THE EST ENTRIES ARE
 
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
 
!     PHIOUT (ESTA) CONTAINS THE FOLLOWING WHERE NGP IS THE NUMBER
!     OF GRID POINTS
 
!     ELEMENT ID
!     NGP SIL NUMBERS
!     NGP VALUES OF THE SHAPE FUNCTIONS AT THIS STRESS POINT
!     REFERENCE TEMPERATURE
!     6 THERMAL STRESS COEFFICIENTS
!     NGP, 6 BY 3 MATRICES, RELATING STRESS TO DISPLACEMENTS AT THIS
!          STRESS POINT (STORED ROW-WISE)
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 INTEGER, INTENT(OUT)                     :: strspt
 INTEGER, INTENT(OUT)                     :: nip
 LOGICAL :: tdep     ,anis       ,rect       ,mtdep
 INTEGER :: cid      ,bgpid      , iest(1)    ,  &
     eid      ,iphio(1)   , itab(3,64) , ib(46)
 REAL :: nu       ,jacob      ,dshpb(3,32),bxyz(3)    ,  &
     gauss(8) ,s(4)       ,gmat(36)   ,store(18)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 COMMON /system/ sysbuf   ,iprnt      ,junk(7)    ,mtemp
 COMMON /matin / mid      ,inflag     ,temp
 COMMON /matout/ e        ,g          ,nu         ,rho        ,  &
     alpha    ,tref       ,SPACE(19)  ,mtdep
 COMMON /matiso/ bufm6(46)
 COMMON /sdr2x5/ est(100) ,phiout(649)
 COMMON /sdr2x6/ cid      ,bgpid(32)  ,eid        ,bgpdt(3,32),  &
     gpt(32)  ,jacob(3,3) ,dshp(3,32) ,detj       ,  &
     d        ,e1         ,e2         ,e3         ,  &
     t(3,3)   ,ngp        ,sglob(18)
 EQUIVALENCE     (est(1),iest(1),dshpb(1,1)),  &
     (phiout(1),iphio(1)),(est(97),idxyz), (est(98),bxyz(1))   ,(ib(1),bufm6(1))
 DATA    gauss/  .57735027, .55555556, .77459667, .88888889 ,  &
     .34785485, .86113631, .65214515, .33998104 /
 
 IF (strspt == 0) strspt = strspt + 1
 IF (strspt > 1) GO TO 505
 
!     MOVE EST DATA INTO /SDR2X6/, /MATIN/, AND PHIOUT
 
 eid = iest(1)
 ngp = 12*TYPE - 4
 mid = iest(ngp+2)
 cid = iest(ngp+3)
 nip = iest(ngp+4)
 IF (nip == 0) nip = TYPE/2 + 2
 
!     FOR STRESS COMPUTATION, SET NUMBER OF STRESS POINTS TO 2
!     NUMBER OF GAUSS POINTS) TO CUT DOWN ON AMOUNT OF INFO ON ESTA
 
 nip = 2
 l   = 0
 DO  i = 1,nip
   DO  j = 1,nip
     DO  k = 1,nip
       l   = l + 1
       itab(1,l) = i
       itab(2,l) = j
       itab(3,l) = k
     END DO
   END DO
 END DO
 DO  i = 1,ngp
   gpt(  i) = est (5*ngp+7+i)
   bgpid(i) = iest(ngp+4+4*i)
   DO  j = 1,3
     bgpdt(j,i) = est(ngp+4+4*i+j)
   END DO
 END DO
 phiout(1) = est(1)
 DO  i = 1,ngp
   phiout(i+1) = est(i+1)
 END DO
 
!     FETCH MATERIAL PROPERTIES
 
!     CHANGE FOR GENERAL ANISOTROPIC MATERIAL
 
!     TEST FOR ANISOTROPIC MATERIAL
 
 anis   = .false.
 inflag = 10
 
!     TEST FOR RECTANGULAR COORDINATE SYSTEM IN WHICH ANISOTROPIC
!     MATERIAL IS DEFINED
 
 rect = .true.
 tdep = .true.
 
 DO  i = 2,ngp
   IF (gpt(i) /= gpt(1)) GO TO 70
 END DO
 tdep = .false.
 70 temp = gpt(1)
 CALL mat (eid)
 IF (ib(46) == 6) anis = .true.
 tref = bufm6(44)
 IF (.NOT.mtdep) tdep = .false.
 
!     IF ISOTROPIC, TEMPERATURE INDEPENDENT MATERIAL, COMPUTE CONSTANTS
 
 IF (tdep) GO TO 500
 IF (anis) GO TO 490
 IF (ib(46) /= 0) GO TO 480
 WRITE  (iprnt,470) uwm,mid,eid
 470 FORMAT (a25,' 4005. AN ILLEGAL VALUE OF -NU- HAS BEEN SPECIFIED',  &
     ' UNDER MATERIAL ID =',i10,' FOR ELEMENT ID =',i10, /32X,  &
     'NU = 0.333 ASSUMED FOR STRESS COMPUTATION')
 e1 = 1.5*e
 e2 = 0.75*e
 e3 = 0.375*e
 GO TO 490
 480 e1 = bufm6(1)
 e2 = bufm6(2)
 e3 = bufm6(22)
 alpha = bufm6(38)
 GO TO 500
 
!     IF MATERIAL IS ANISOTROPIC, DEFINED IN A RECTANGULAR
!     COORDINATE SYSTEM, AND NOT TEMPERATURE DEPENDENT, TRANSFORM
!     IT TO THE BASIC SYSTEM.
 
 490 IF (.NOT.rect) GO TO 500
 
!     ADD CODE TO TRANSFORM GENERAL ANISOTROPIC MATERIAL
!     TO BASIC COORDINATE SYSTEM HERE.
 
 DO  ijk = 1,36
   gmat(ijk) = bufm6(ijk)
 END DO
 
!     INITIALIZATION TO FIND GAUSS POINT COORDINATES
 
 505 CONTINUE
 500 nipm1 = nip - 1
 SELECT CASE ( nipm1 )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
 END SELECT
 510 s(1) = gauss(1)
 s(2) =-gauss(1)
 GO TO  540
 520 s(1) = gauss(3)
 s(2) = 0.
 s(3) =-gauss(3)
 GO TO  540
 530 s(1) = gauss(6)
 s(2) = gauss(8)
 s(3) =-gauss(8)
 s(4) =-gauss(6)
 540 IF (strspt == nip**3+1) GO TO 541
 l = itab(1,strspt)
 x = s(l)
 l = itab(2,strspt)
 y = s(l)
 l = itab(3,strspt)
 z = s(l)
 GO TO 542
 541 x = 0.
 y = 0.
 z = 0.
 542 CONTINUE
 
!     GENERATE SHAPE FUNCTIONS AND JACOBIAN MATRIX INVERSE
 
 CALL ihexss (TYPE,phiout(ngp+2),dshp,jacob,detj,eid,x,y,z,bgpdt)
 IF (detj /= 0.0) GO TO 605
 
!     FALL HERE IF JACOBIAN MATRIX SINGULAR (BAD ELEMENT)
 
 j = ngp*19 + 7
 DO  i = 1,j
   phiout(ngp+1+i) = 0.0
 END DO
 RETURN
 
!     COMPUTE STRAIN-DISPLACEMENT RELATIONS
 
!     REVERSE CALLING SEQUENCE SINCE MATRICES ARE COLUMN STORED
 
 605 CALL gmmats (dshp,ngp,3,0,jacob,3,3,0,dshpb)
 
!     IF MATERIAL IS TEMPERATURE DEPENDENT, MUST COMPUTE TEMPERATURE
!     AT THIS STRESS POINT AND FETCH MATERIAL PROPERTIES AGAIN
 
 IF (.NOT.tdep) GO TO 620
 temp = 0.0
 DO  j = 1,ngp
   temp = temp + gpt(j)*phiout(ngp+1+j)
 END DO
 CALL mat (eid)
 IF (anis) GO TO 630
 IF (ib(46) /= 0) GO TO 615
 WRITE (iprnt,470) uwm,mid,eid
 e1 = 1.5*e
 e2 = 0.75*e
 e3 = 0.375*e
 GO TO 640
 615 e1 = bufm6(1)
 e2 = bufm6(2)
 e3 = bufm6(22)
 alpha = bufm6(38)
 GO TO 640
 
!     IF MATERIAL IS ANISOTROPIC AND NOT DEFINED IN RECTANGJLAR
!     COORDINATE SYSTEM, TRANSFORM IT TO BASIC COORDINATE SYSTEM AT
!     THIS STRESS POINT.
 
 
!     IN THIS VERSION, ANISOTROPIC PROPERTIES MUST BE RECTANGULAR
!     JUST STORE G MATRIX
!     ===========================================================
 
!     THIS CODE MUST BE COMPLETED WHEN GENERAL ANISOTROPIC MATERIAL IS
!     ADDED.
 
 620 IF (.NOT.anis) GO TO 640
 630 CONTINUE
 DO  ijk = 1,36
   gmat(ijk) = bufm6(ijk)
 END DO
 
!     INSERT GLOBAL TO BASIC TRANSFORMATION OPERATIONS HERE FOR
!     ANISOTROPIC MATERIAL.
 
!     MATERIAL HAS BEEN EVALUATED AT THIS STRESS POINT WHEN GET TO HERE
 
!     TEMPERATURE TO STRESS VECTOR
 
 640 phiout(2*ngp+2) = tref
 IF (anis) GO TO 660
 
!     ISOTROPIC CASE
 
 DO  j = 1,3
   phiout(2*ngp+2+j) = -alpha*(e1+2.0*e2)
   phiout(2*ngp+5+j) = 0.0
 END DO
 GO TO 670
 
!     ANISOTROPIC CASE
 
!     ADD CODE WHEN ANISOTROPIC MATERIAL BECOMES AVAILABLE
 
 660 CONTINUE
 CALL gmmats (gmat,6,6,0,bufm6(38),6,1,0,phiout(2*ngp+3))
 DO  ijk = 1,6
   is = 2*ngp + 2 + ijk
   phiout(is) = -phiout(is)
 END DO
 
!     DISPLACEMENT TO STRESS MATRICES
 
 670 DO  i = 1,ngp
   is = 2*ngp + 8 + 18*(i-1)
   
!     ROW-STORED
   
   IF (anis) GO TO 680
   
!     ISOTROPIC CASE
   
   phiout(is+ 1) = e1*dshpb(1,i)
   phiout(is+ 2) = e2*dshpb(2,i)
   phiout(is+ 3) = e2*dshpb(3,i)
   phiout(is+ 4) = e2*dshpb(1,i)
   phiout(is+ 5) = e1*dshpb(2,i)
   phiout(is+ 6) = e2*dshpb(3,i)
   phiout(is+ 7) = e2*dshpb(1,i)
   phiout(is+ 8) = e2*dshpb(2,i)
   phiout(is+ 9) = e1*dshpb(3,i)
   phiout(is+10) = e3*dshpb(2,i)
   phiout(is+11) = e3*dshpb(1,i)
   phiout(is+14) = e3*dshpb(3,i)
   phiout(is+15) = e3*dshpb(2,i)
   phiout(is+16) = e3*dshpb(3,i)
   phiout(is+18) = e3*dshpb(1,i)
   phiout(is+12) = 0.0
   phiout(is+13) = 0.0
   phiout(is+17) = 0.0
   GO TO 690
   
!     ANISOTROPIC CASE
   
!     ADD CODE WHEN GENERAL ANISOTROPIC MATERIAL BECOMES AVAILABLE
   
   680 CONTINUE
   DO  ijk = 1,18
     store(ijk) = 0.
   END DO
   store( 1) = dshpb(1,i)
   store( 5) = dshpb(2,i)
   store( 9) = dshpb(3,i)
   store(10) = dshpb(2,i)
   store(11) = dshpb(1,i)
   store(14) = dshpb(3,i)
   store(15) = dshpb(2,i)
   store(16) = dshpb(3,i)
   store(18) = dshpb(1,i)
   
   CALL gmmats (gmat(1),6,6,0,store(1),6,3,0,phiout(is+1))
   
!     POST-MULTIPLY BY GLOBAL TO BASIC TRANSFORMATION MATRIX,
!     IF NECESSARY
   
   690 IF (bgpid(i) == 0) CYCLE
   idxyz = bgpid(i)
   DO  k = 1,3
     bxyz(k) = bgpdt(k,i)
   END DO
   
!     FETCH TRANSFORMATION AND USE IT
   
   CALL transs (idxyz,t)
   CALL gmmats (phiout(is+1),6,3,0,t,3,3,0,sglob)
   DO  j = 1,18
     phiout(is+j) = sglob(j)
   END DO
 END DO
 iphio(20*ngp+9) = nip
 nwdnow = 20*ngp + 9
 nwdiso = 649 - nwdnow
 IF (nwdiso == 0) RETURN
 DO  i = 1,nwdiso
   isub = nwdnow + i
   phiout(isub) = 0.
 END DO
 RETURN
END SUBROUTINE sihex1
