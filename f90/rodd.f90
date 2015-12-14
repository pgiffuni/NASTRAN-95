SUBROUTINE rodd
     
!     THIS ROUTINE PROCESSES ROD ELEMENT DATA TO PRODUCE STIFFNESS AND
!     MASS MATRICES. IF THE HEAT TRANSFER OPTION IS ON, CONDUCTIVITY AND
!     CAPACITY MATRICES ARE PRODUCED
 
!     THIS ROUTINE CAN COMPUTE BOTH CONVENTIONAL AND CONSISTENT
!     MASS MATRICES
 
!     DOUBLE PRECISION VERSION
 
!     THIS VERSION WAS SPECIALLY CODED TO ILLUSTRATE A GENERAL
!     USE OF THE IMPROVED MATRIX GENERATOR.
 
!     THE EST ENTRY FOR THIS ELEMENT CONTAINS
 
!     POSITION     NAME       DESCRIPTION
!     *****        *****      *******************************
!     1             EID       ELEMENT ID NO.
!     2             SIL1      SCALAR INDEX OF POINT A
!     3             SIL2      SCALAR INDEX OF POINT B
!     4             MID       MATERIAL DATA ID
!     5             AFACT     AREA OF CROSS SECTION
!     6             JFACT     TORSIONAL STIFFNESS COEFFICIENT
!     7             CFACT     TORSIONAL STRESS RECOVERY DISTANCE
!     8             MU        NON-STRUCTURAL MASS PER LENGTH
!     9-16          BGPDT     BASIC GRID POINT DATA. COORDINATE SYSTEM
!                             NUMBER AND  X,Y,Z LOCATION FOR 2 POINTS
!     17            TBAR      AVERAGE ELEMENT TEMPERATURE
 
 
 LOGICAL :: nogo
 INTEGER :: sil1     ,sil2     ,iest(13) ,eid      ,GE       ,  &
     dict(7)  ,elid     ,estid
 REAL :: jfact    ,mu       ,kcon     ,est(200)
 DOUBLE PRECISION :: evect(3) ,el       ,ke       ,me       ,  &
     te       ,ha(3)    ,hb(3)    ,kha(3)   ,khb(3)   ,  &
     ta(9)    ,tb(9)    ,scale    ,k        ,mjidum(9),  &
     massii(9),massjj(9),massij(9),massji(9),mijdum(9)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /matin / matid    ,inflag   ,eltemp   ,dum(3)
 COMMON /matout/ e        ,g        ,nu       ,rho      ,alfa     ,  &
     tsub0    ,GE
 COMMON /hmtout/ kcon
 COMMON /emgprm/ ixtra    ,izr      ,nzr      ,dumy(12) ,kmbgg(3) ,  &
     iprec    ,nogo     ,heat     ,icmbar
 COMMON /emgdic/ dum2(2)  ,nlocs    ,elid     ,estid
 COMMON /zzzzzz/ k(1)
 
!     THE VARIABLE K IS OPEN CORE. OPEN SPACE EXISTS FROM Z(IZ) TO Z(NZ)
!     THIS IS INTENDED AS AN EXAMPLE. NORMALLY FOR SMALL ARRAYS
!     LOCAL VARIABLES MAY BE USED.
 
 COMMON /emgest/ eid      ,sil1     ,sil2     ,mid      ,afact    ,  &
     jfact    ,cfact    ,mu       ,bgpdt(4,2),tbar
 COMMON /system/ ksystm(63)
 EQUIVALENCE     (ksystm( 2),ioutpt),(ksystm(56),iheat)  ,  &
     (eid,est(1),iest(1)),(cp,kcon)
 
!     FOR DOUBLE PRECISION THE POINTERS TO OPEN CORE MUST BE MODIFIED.
 
 iz = (izr-2)/iprec + 2
 nz = nzr/iprec
 IF (nz-iz <= 144) GO TO 290
 dict(1) = estid
 
!     SUBTRACT BASIC LOCATIONS TO OBTAIN LENGTH ETC.
 
 DO  i = 1,3
   evect(i) = bgpdt(i+1,2) - bgpdt(i+1,1)
 END DO
 
 el = DSQRT(evect(1)**2 + evect(2)**2 + evect(3)**2)
 IF (el <= 0.0D0) GO TO 270
 
!     IF HEAT TRANSFER PROBLEM TRANSFER.  CALL MATERIAL SUBROUTINE
 
 inflag = 1
 matid  = mid
 eltemp = tbar
 IF (iheat == 1) GO TO 240
 CALL mat (eid)
 ke = DBLE(e*afact)/el
 me = (DBLE(rho*afact+mu))*el/2.0D0
 te = DBLE(g*jfact)/el
 
!     PROCESS STIFFNESS HERE
 
 IF (kmbgg(1) == 0) GO TO 220
 IF (ke == 0.0D0 .AND. te == 0.0D0) GO TO 220
 
!     GENERATE   HA  =  (E*TA)/EL   AND  HB = (E*TB)/EL
 
 IF (iest(9) == 0) GO TO 30
 CALL transd (bgpdt(1,1),ta)
 CALL gmmatd (evect,1,3,0, ta,3,3,0, ha)
 DO  i = 1,3
   ha(i) = ha(i)/el
 END DO
 GO TO 50
 30 DO  i = 1,3
   ha(i) = evect(i)/el
 END DO
 50 IF (iest(13) == 0) GO TO 70
 CALL transd (bgpdt(1,2),tb)
 CALL gmmatd (evect,1,3,0, tb,3,3,0, hb)
 DO  i = 1,3
   hb(i) = hb(i)/el
 END DO
 GO TO 90
 70 DO  i = 1,3
   hb(i) = evect(i)/el
 END DO
 
!     THE GENERAL 12X12  MATRIX FOR THE ROD ELEMENT IS
!                            -                              -
!                            1HA K HA1   0  1HA K HB1       1
!                **   **     1 ------1------1-------1-------1
!                *  K  *   = 1   0   1HA T A1       1HA T HB1
!                **   **     1 ------1------1-------1-------1
!                            1HB K HA1      1HB K HB1       1
!                            1 ------1------1-------1-------1
!                            1       1HB T A1       1HB T HB1
!                            1       1      1       1       1
!                            -                              -
!                      EACH BLOCK  ABOVE IS A 3 BY 3 MATRIX
 
!     TEST AND SET COMPONENT CODE    111= 7     111000=56
 
 90 icode = 0
 ndof  = 0
 IF (te  /= 0.d0) GO TO 100
 icode = 7
 ndof  = 6
 GO TO 120
 100 IF (ke  /= 0.d0) GO TO 110
 icode = 56
 ndof  = 6
 GO TO 120
 110 icode = 63
 ndof  = 12
 120 nsq   = ndof**2
 ng    = ndof/2
 npart = ng*ndof
 izero = iz - 1
 ipass = 1
 DO   i = 1,nsq
   izpi  = iz + i - 1
   k(izpi) = 0.0D0
 END DO
 
!     EXTENSIONAL STIFFNESS TERMS ARE COMPUTED HERE.
 
 IF (icode == 56) GO TO 200
 scale = ke
 140 DO  i = 1,3
   kha(i) = scale*ha(i)
   khb(i) = scale*hb(i)
 END DO
 
!     THE MATRIX COLUMNS AND ROWS MUST BE IN THE NUMERICAL ORDER
!     OF TH SIL VALUES. THE POINTERS INTO THE MATRIX ARE VARIABLES.
 
 IF (sil2-sil1 < 0.0) THEN
   GO TO   160
 ELSE IF (sil2-sil1 == 0.0) THEN
   GO TO   270
 ELSE
   GO TO   170
 END IF
 160 ibbz = izero
 iabz = izero + ng
 ibaz = izero + npart
 iaaz = ibaz  + ng
 GO TO 180
 170 iaaz = izero
 ibaz = izero + ng
 iabz = izero + npart
 ibbz = iabz  + ng
 180 CONTINUE
 DO  j = 1,3
   DO  i = 1,3
     ij  = ndof*(j-1) + i
     iaa = ij + iaaz
     k(iaa) = kha(i)*ha(j)
     iba = ij + ibaz
     k(iba) =-khb(i)*ha(j)
     iab = ij + iabz
     k(iab) =-kha(i)*hb(j)
     ibb = ij + ibbz
     k(ibb) = khb(i)*hb(j)
   END DO
 END DO
 
!     THE TORSIONAL STIFFNESS TERMS ARE FORMED USING TE INSTEAD OF KE
!     THEY ARE INSERTED IN THE MATRIX WITH  A CONSTANT OFFSET, 3*12+3.
 
 200 IF (ipass == 2) GO TO 210
 IF (ndof == 12) izero = 38 + iz
 ipass = 2
 scale = te
 IF (icode /= 7) GO TO 140
 210 ipart = iz
 dict(2) = 1
 dict(3) = ndof
 dict(4) = icode
 dict(5) = GE
 CALL emgout (k(ipart),k(ipart),nsq,1,dict,1,iprec)
 
!     THE MASS MATRIX TERMS ARE CALCULATED HERE.
 
 220 IF (kmbgg(2) == 0 .OR. me == 0.0D0) RETURN
 dict(3) = 6
 dict(4) = 7
 dict(5) = 0
 
!     CHECK TO SEE IF CONVENTIONAL OR CONSISTENT MASS MATRIX IS REQUIRED
 
 IF (icmbar > 0) GO TO 400
 
!     CONVENTIONAL MASS MATRIX TERMS ARE COMPUTED HERE
 
 dict(2) = 2
 ldata   = 6
 izp5    = iz + 5
 DO  i = iz,izp5
   k(i) = me
 END DO
 GO TO 600
 
!     CONSISTENT MASS MATRIX TERMS ARE COMPUTED HERE
 
 400 dict(2) = 1
 ldata   = 36
 DO  i = 1,9
   massii(i) = 0.0D0
   massjj(i) = 0.0D0
   massij(i) = 0.0D0
   massji(i) = 0.0D0
   mijdum(i) = 0.0D0
   mjidum(i) = 0.0D0
 END DO
 me = 2.0D0*me
 DO  i = 1,9,4
   massii(i) = me/3.0D0
   massjj(i) = me/3.0D0
   massij(i) = me/6.0D0
   massji(i) = me/6.0D0
   mijdum(i) = me/6.0D0
   mjidum(i) = me/6.0D0
 END DO
 IF (sil2-sil1 < 0.0) THEN
   GO TO   480
 ELSE IF (sil2-sil1 == 0.0) THEN
   GO TO   270
 END IF
 460 iti = 9
 itj = 13
 GO TO 500
 480 iti = 13
 itj = 9
 500 IF (iest(iti) == 0) GO TO 520
 CALL transd (iest(iti), ta)
 CALL gmmatd (ta,3,3,1, massii,3,3,0, k(iz))
 CALL gmmatd (k(iz),3,3,0, ta,3,3,0, massii)
 CALL gmmatd (ta,3,3,1, mijdum,3,3,0, massij)
 CALL gmmatd (mjidum,3,3,0, ta,3,3,0, massji)
 520 IF (iest(itj) == 0) GO TO 560
 CALL transd (iest(itj), ta)
 CALL gmmatd (ta,3,3,1, massjj,3,3,0, k(iz))
 CALL gmmatd (k(iz),3,3,0, ta,3,3,0, massjj)
 CALL gmmatd (massij,3,3,0, ta,3,3,0, mijdum)
 CALL gmmatd (ta,3,3,1, massji,3,3,0, mjidum)
 DO  i = 1,9
   massij(i) = mijdum(i)
   massji(i) = mjidum(i)
 END DO
 560 DO  i = 1,3
   kz = iz + i - 1
   k(kz   ) = massii(i  )
   k(kz+ 6) = massii(i+3)
   k(kz+12) = massii(i+6)
   k(kz+ 3) = massij(i  )
   k(kz+ 9) = massij(i+3)
   k(kz+15) = massij(i+6)
   k(kz+18) = massji(i  )
   k(kz+24) = massji(i+3)
   k(kz+30) = massji(i+6)
   k(kz+21) = massjj(i  )
   k(kz+27) = massjj(i+3)
   k(kz+33) = massjj(i+6)
 END DO
 600 CALL emgout (k(iz),k(iz),ldata,1,dict,2,iprec)
 RETURN
 
!     HEAT TRANSFER CALCULATIONS ARE PERFORMED HERE
 
 240 inflag  = 1
 dict(2) = 1
 dict(3) = 2
 dict(4) = 1
 dict(5) = 0
 IF (kmbgg(1) == 0) GO TO 250
 CALL hmat (eid)
 k(iz) = DBLE(afact*kcon)/el
 IF (k(iz) == 0.0D0) GO TO 250
 k(iz+1) = -k(iz)
 k(iz+2) = -k(iz)
 k(iz+3) =  k(iz)
 CALL emgout (k(iz),k(iz),4,1,dict,1,iprec)
 250 inflag = 4
 IF (kmbgg(1) == 0) RETURN
 CALL hmat (eid)
 k(iz) = DBLE(afact*cp)*el/2.0D0
 IF (k(iz) == 0.0D0) RETURN
 k(iz+1) = k(iz)
 dict(2) = 2
 CALL emgout (k(iz),k(iz),2,1,dict,3,iprec)
 RETURN
 
 270 nogo = .true.
 WRITE  (ioutpt,280) ufm,eid
 280 FORMAT (a23,' 3118, ROD ELEMENT NO.',i9,  &
     ' HAS ILLEGAL GEOMETRY OR CONNECTIONS.')
 RETURN
 
 290 nogo = .true.
 WRITE  (ioutpt,300) ufm
 300 FORMAT (a23,' 3119, INSUFFICIENT CORE TO PROCESS ROD ELEMENTS')
 RETURN
END SUBROUTINE rodd
