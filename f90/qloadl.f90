SUBROUTINE qloadl (iopt)
     
!     THIS ROUTINE CALCULATES THERMAL LOADS FROM QBDY1, QBDY2, OR
!     QVECT DATA. THE INPUT DATA, READ FROM FILE SLT, IS -
 
!     ENTRY       QBDY1         QBDY2          QVECT
!     -----       -----         -----          -----
!       1          TYPE         EL.ID.          SIL1
!       2         EL.ID.         TYPE           SIL2
!       3          SIL1          SIL1           SIL3
!       4          SIL2          SIL2           SIL4
!       5          SIL3          SIL3          EL.ID.
!       6          SIL4          SIL4           TYPE
!      7-10      C1,C2,C3,C4    -SAME          -SAME
!     11-13       E VECTOR       NONE           NONE
!     14-16      V1 VECTOR         *             *
!     17-19      V2 VECTOR         *             *
 
 
 INTEGER, INTENT(IN OUT)                  :: iopt
 LOGICAL :: nogo,transt
 INTEGER :: slt,OLD,bg,subr(2),igrids(6),TYPE,sils(4),ie(3), minus(2)
 REAL :: e(3),coef(4),v1(3),v2(3),card(19)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /loadx / lc,slt,bg,OLD,nxn(12),ifm,nyn(2),ilid
 COMMON /zzzzzz/ core(1)
 COMMON /system/ sysbuf,iout
 COMMON /qvect / itran,iqvect
 EQUIVALENCE     (TYPE,card(1)),(id,card(2)),(sils(1),card(3)),  &
     (coef(1),card(7)),(e(1),ie(1),card(11)), (v1(1),card(14)),(v2(1),card(17))
 DATA    igrids/ 1,2,2,3,4,2   /
 DATA    subr  / 4HQLOA,4HDL   /
 DATA    itran1, iold,minus    / 4HTRAN,0,-1,-1 /
 
 transt = .false.
 IF (itran == itran1) transt = .true.
 nwords = 10
 IF (iopt == 3) nwords = 19
 
 CALL READ (*100,*110,slt,card(1),nwords,0,flag)
 
!     REARRANGE CARD ARRAY FOR UNIFORMITY.
 
 SELECT CASE ( iopt )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 10
   CASE (    3)
     GO TO 40
 END SELECT
 10 dot = card(1)
 card(1) = card(2)
 card(2) = dot
 20 n = igrids(TYPE)
 
!     QBDY1 OR QBDY2
 
 DO  i = 1,n
   isil = sils(i)
   core(isil) = core(isil) + coef(i)
 END DO
 RETURN
 
!     QVECT LOADS
 
 40 dot     = card(5)
 dot2    = card(6)
 card(6) = card(4)
 card(5) = card(3)
 card(4) = card(2)
 card(3) = card(1)
 card(2) = dot
 card(1) = dot2
 n       = igrids(TYPE)
 dot     = 0.0
 INT     = 0
 IF (TYPE == 6) GO TO 70
 DO  i = 1,3
   IF (numtyp(ie(i)) == 1) GO TO 51
   dot = dot + e(i)*v1(i)
   CYCLE
   51 INT = INT + 1
 END DO
 IF (INT >   0) GO TO 90
 IF (dot >= 0.0) RETURN
 DO  i = 1,n
   isil = sils(i)
   core(isil) = core(isil) - dot*coef(i)
 END DO
 RETURN
 
!     QVECT ON ELCYL ELEMENT
 
 70 dot2 = 0.0
 DO  i = 1,3
   IF (numtyp(ie(i)) == 1) GO TO 81
   dot  = dot  + e(i)*v1(i)
   dot2 = dot2 + e(i)*v2(i)
   CYCLE
   81 INT = INT + 1
 END DO
 IF (INT > 0) GO TO 90
 coef(1) = coef(1)*SQRT(dot**2 + dot2**2)
 coef(2) = coef(1)
 isil = sils(1)
 core(isil) = core(isil) + coef(1)
 isil = sils(2)
 core(isil) = core(isil) + coef(2)
 RETURN
 
!     GOES HERE IF INTEGERS ARE FOUND IN E VECTOR
 
 90 IF (.NOT. transt) GO TO 120
 
!     BUILD QVECT RECORDS FOR TRANSIENT
 
 IF (ilid == iold) GO TO 91
 IF (iold ==    0) GO TO 92
 
!     TERMINATE OLD RECORD
 
 CALL WRITE (iqvect,minus,2,0)
 92 iold = ilid
 CALL WRITE (iqvect,ilid,1,0)
 
!     DUMP DATA ON IQVECT
 
 91 CALL WRITE (iqvect,n,1,0)
 DO  i = 1,n
   CALL WRITE (iqvect,sils(i),1,0)
   CALL WRITE (iqvect,coef(i),1,0)
 END DO
 CALL WRITE (iqvect,ie,3,0)
 CALL WRITE (iqvect,v1,6,0)
 RETURN
 
 100 CALL mesage (-2,slt,subr)
 110 CALL mesage (-3,slt,subr)
 120 nogo = .true.
 WRITE  (iout,130) ufm,id
 130 FORMAT (a23,' 3080, ERROR IN QVECT DATA, INTEGER VALUES SPECIFIED'  &
     ,      ' FOR THERMAL FLUX VECTOR COMPONENTS', /30X,  &
     'IN A NON-TRANSIENT ANALYSIS.', /30X,'ELEMENT ID = ',i9)
 CALL mesage (-61,0,subr)
 RETURN
END SUBROUTINE qloadl
