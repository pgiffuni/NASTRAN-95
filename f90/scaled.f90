SUBROUTINE scaled (TYPE,emord)
     
!     THIS ROUTINE PROCESSES CELAS, CDAMP, AND CMASS ELEMENTS.
 
!     TYPE  - DENOTES FORM OF EST DATA. IE CELAS1,CELAS2,ETC.
!     EMORD - DENOTES MATRIX  1 = CELAS = STIFFNESS MATRIX,
!                             2 = CMASS = MASS MATRIX,
!                             3 = CDAMP = DAMPING MATRIX
 
!     EST FOR ELAS ELEMENTS
 
!                     TYPE           TYPE           TYPE           TYPE
!             CELAS1         CELAS2         CELAS3         CELAS4
!             ------  ----   ------  ----   ------  ----   ------  ----
!     ECPT(1) IELID     I    IELID     I    IELID     I    IELID     I
!     ECPT(2) IGP1      I    K         R    IS1       I    K         R
!     ECPT(3) IGP2      I    IGP1      I    IS2       I    IS1       I
!     ECPT(4) IC1       I    IGP2      I    K         R    IS2       I
!     ECPT(5) IC2       I    IC1       I    GSUBE     R
!     ECPT(6) K         R    IC2       I    S         R
!     ECPT(7) GSUBE     R    GSUBE     R
!     ECPT(8) S         R    S         R
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 INTEGER, INTENT(IN OUT)                  :: emord
 LOGICAL :: nogo
 INTEGER :: eid,isil(2),icomp(2),gpt(4),cpt(2),  &
     kpt(4),gspt(4),code,iest(1),dict(7),gsube,elid, estid
 DOUBLE PRECISION :: dz(16)
 DIMENSION        z(16)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm
 COMMON /emgest/  est(100)
 COMMON /emgprm/  dumy(15),imat(3),iprec,nogo
 COMMON /system/  ksystm(65)
 COMMON /emgdic/  dum2(2),nlocs,elid,estid
 EQUIVALENCE      (ksystm(2),ioutpt),(z(1),dz(1)),(iest(1),est(1))
 DATA    gpt   /  2, 3, 2, 3 /, cpt / 4, 5 /, kpt /6, 2, 4, 2 /
 DATA    gspt  /  7, 7, 5, 0 /
 
!     TEST IF MATRIX TO BE PRODUCED IS REQUESTED
 
 IF (imat(emord) == 0) RETURN
 
!     MOVE EST DATA TO LOCAL ARRAYS.  LOCATIONS ARE GIVEN BY DATA //
 
 eid     = iest(1)
 ip      = kpt(TYPE)
 z(1)    = est(ip)
 gsube   = 0
 icomp(1)= 0
 icomp(2)= 0
 dict(2) = 1
 ngrids  = 2
 ip      = gpt(TYPE)
 isil(1) = iest(ip)
 isil(2) = iest(ip+1)
 IF (TYPE >= 3) GO TO 10
 ip = cpt(TYPE)
 IF (iest(ip  ) /= 0) icomp(1) = iest(ip  ) - 1
 IF (iest(ip+1) /= 0) icomp(2) = iest(ip+1) - 1
 
!     IF ONE SIL IS ZERO INSURE THAT IT IS THE SECOND.
!     IF BOTH SILS ARE NON-ZERO MAKE SURE HIGHER OF TWO IS SECOND.
 
 10 IF (isil(2) == 0) GO TO 5
 IF (isil(1) == 0) GO TO 4
 IF (isil(1) <= isil(2)) GO TO 5
 
!     SWITCH SILS AND COMPS
 
 4 ip      = isil(2)
 isil(2) = isil(1)
 isil(1) = ip
 ip      = icomp(2)
 icomp(2)= icomp(1)
 icomp(1)= ip
 5 IF (isil(2) > 0) GO TO 20
 
!     IF THE SECOND SIL EQUALS ZERO THE ELEMENT IS GROUNDED
!     ONLY A SINGLE MATRIX TERM IS PRODUCED
 
 ngrids = 1
 dict(2)= 1
 nterms = 1
 code   = 2**icomp(1)
 ncol   = 1
 GO TO 80
 
 20 IF (isil(2) /= isil(1)) GO TO 30
 
!     IF THE ELEMENT CONNECTS TWO COMPONENTS OF THE SAME POINT IT
!     MUST HAVE SPECIAL TREATMENT
 
 IF (icomp(2) == icomp(1)) GO TO 110
 
!     IN THE GENERAL CASE, THE CONNECTED COMPONENTS MAY BE THE SAME
!     AND THE MATRIX IS A 2 BY 2.  IF THE COMPONENTS ARE DIFFERENT
!     THE MATRIX WILL BE A 4 BY 4 WITH ADDITIONAL ZEROS.
 
 GO TO 40
 30 IF (icomp(1) == icomp(2)) GO TO 70
 
 40 nterms= 16
 code  = 2**icomp(1) + 2**icomp(2)
 ncol  = 4
 DO  i = 2,16
   z( i) = 0.0
 END DO
 IF (icomp(2) < icomp(1)) GO TO 60
 z( 4) =-z(1)
 z(13) =-z(1)
 z(16) = z(1)
 IF (isil(1) /= isil(2)) GO TO 80
 z( 2) = z( 4)
 z( 5) = z(13)
 z( 6) = z(16)
 z( 4) = 0.0
 z(13) = 0.0
 z(16) = 0.0
 GO TO 80
 60 z( 6) = z(1)
 z( 7) =-z(1)
 z(10) =-z(1)
 z(11) = z(1)
 z( 1) = 0.0
 IF (isil(1) /= isil(2)) GO TO 80
 z( 1) = z(11)
 z( 2) = z(10)
 z( 5) = z( 7)
 z( 7) = 0.0
 z(10) = 0.0
 z(11) = 0.0
 GO TO 80
 
!     COMPONENTS ARE THE SAME FOR BOTH POINTS
 
 70 nterms= 4
 ncol  = 2
 code  = 2**icomp(1)
 z(2)  =-z(1)
 z(3)  =-z(1)
 z(4)  = z(1)
 
!     OUTPUT THE MATRIX HERE
 
 80 dict(1) = estid
 dict(3) = ncol
 dict(4) = code
 dict(5) = 0
 ipg     = gspt(TYPE)
 
!     STRUCTURAL DAMPING FOR  STIIFNESS MATRICES IS INSERTED IN DICT
 
 IF (emord == 1 .AND. TYPE <= 3) dict(5) = iest(ipg)
 IF (iprec == 1) GO TO 100
 i = nterms
 90 dz(i) = z(i)
 i = i - 1
 IF (i > 0) GO TO 90
 100 CALL emgout (z,dz,nterms,1,dict,emord,iprec)
 RETURN
 
 110 WRITE  (ioutpt,120) uwm,eid
 120 FORMAT (a25,' 3120, IMPROPER CONNECTION ON CELAS ELEMENT',i9)
 RETURN
END SUBROUTINE scaled
