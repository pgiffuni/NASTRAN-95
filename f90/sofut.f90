SUBROUTINE sofut
     
!     THE PURPOSE OF THE MODULE IS TO PERFORM THE TASKS OF ALTERING THE
!     SOF FILE IN ORDER TO EDIT, PURGE, AND EQUIVALENCE THE DATA ITEMS
!     OF SELECTED SUBSTRUCTURES.  THE CALLING SEQUENCE TO THE MODULE IS
 
!     SOFUT     //V,N,DRY/C,N,NAME1/C,N,OPER/C,N,OPT/C,N,NAME2/
!                 C,N,PREFX/C,N,IA/C,N,IB/C,N,IC/C,N,ID/C,N,IE $
 
 EXTERNAL        rename
 LOGICAL :: ditup
 INTEGER :: dry,oper,opt,prefx,sysbuf,dele,renam,NAME(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / dry,name1(2),oper(2),opt,name2(2),prefx(2), items(10)
 COMMON /sof   / sss(33),ditup
 COMMON /zzzzzz/ iz(1)
 COMMON /system/ sysbuf,nout
 DATA    iedit , idest,iequiv / 4HEDIT    ,4HDEST ,4HEQUI   /
 DATA    iprnt / 4HSOFP/
 DATA    dele  / 4HDELE/
 DATA    renam / 4HRENA/
 DATA    NAME  / 4HSOFU,4HT   /
 DATA    iscr1 / 301   /
 
 itask = 0
 IF (oper(1) ==  iedit) itask = 1
 IF (oper(1) ==  idest) itask = 2
 IF (oper(1) == iequiv) itask = 3
 IF (oper(1) ==  iprnt) itask = 4
 IF (oper(1) ==   dele) itask = 5
 IF (oper(1) ==  renam) itask = 6
 IF (itask == 0) GO TO 1000
 
!     ALLOCATE BUFFERS FOR THE SOF UTILITY SUBROUTINES
 
 nz  = korsz(iz)
 IF (3*sysbuf > nz) CALL mesage (-8,0,NAME(1))
 ib1 = nz  - sysbuf + 1
 ib2 = ib1 - sysbuf - 1
 ib3 = ib2 - sysbuf
 CALL sofopn (iz(ib1),iz(ib2),iz(ib3))
 nz  = ib3 - 1
 SELECT CASE ( itask )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 130
   CASE (    5)
     GO TO 180
   CASE (    6)
     GO TO 200
 END SELECT
 
!     EDIT OPERATION
 
 20 CALL edit (name1(1),opt,itest)
 GO TO 50
 
!     DESTROY OPERATION
 
 30 i = nz/2 + 1
 CALL dstroy (name1(1),itest,iz,iz(i),i-1)
 GO TO 50
 
!     EQUIVALENCE OPERATION
 
 40 i = nz/2 + 1
 CALL seteq (name1,name2,prefx,dry,itest,iz,i-1)
 
!     TEST RETURN CODE
 
 50 SELECT CASE ( itest )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 110
   CASE (    4)
     GO TO 60
   CASE (    5)
     GO TO 110
   CASE (    6)
     GO TO 70
   CASE (    7)
     GO TO 110
   CASE (    8)
     GO TO 80
   CASE (    9)
     GO TO 90
   CASE (   10)
     GO TO 100
 END SELECT
 60 WRITE (nout,1010) uwm,name1
 GO TO 100
 70 WRITE (nout,1020) uwm,name1
 GO TO 100
 80 WRITE (nout,1030) uwm,name2
 GO TO 100
 90 WRITE (nout,1040) uwm,name2
 100 dry = -2
 110 CALL sofcls
 GO TO 1100
 
!     PRINT OPERATIONS
 
 130 IF (opt > 0.0) THEN
   GO TO   150
 END IF
 
!     PRINT SOF TABLE OF CONTENTS (DIT MDI)
 
 140 CALL softoc
 IF (opt == 0) GO TO 170
 
!     PRINT SOF DATA ITEMS
 
 150 DO  i = 1,5
   ii = ittype(items(2*i-1))
   IF (ii < 0) THEN
     GO TO   160
   ELSE IF (ii == 0) THEN
     GO TO   152
   ELSE
     GO TO   154
   END IF
   
!     TABLE ITEM
   
   152 CALL itmprt (name1,items(2*i-1),nz,opt)
   CYCLE
   
!     MATRIX ITEM
   
   154 CALL matwrt (iscr1,name1,items(2*i-1),nz)
   160 CONTINUE
 END DO
 170 CALL sofcls
 GO TO 1100
 
!     DELETE OPERATION
 
 180 DO  i = 1,10
   CALL DELETE (name1,items(i),itest)
 END DO
 GO TO 50
 
!     RENAME OPERATION
 
 200 CALL rename (name1,name2,iz(1),nz,itest)
 GO TO 50
 
!     ERROR MESSAGES
 
 1000 WRITE  (nout,1001) uwm,oper(1),oper(2)
 1001 FORMAT (a25,' 6217, MODULE SOFUT - ',2A4,' IS AN ILLEGAL ',  &
     'PARAMETER NAME.')
 GO TO 1100
 
 1010 FORMAT (a25,' 6212, MODULE SOFUT - THE SUBSTRUCTURE ',2A4,  &
     ' DOES NOT EXIST.')
 
 1020 FORMAT (a25,' 6218, MODULE SOFUT - THE SUBSTRUCTURE ',2A4,1X,  &
     'CANNOT BE DESTROYED BECAUSE IT IS AN IMAGE SUBSTRUCTURE.')
 
 1030 FORMAT (a25,' 6219, MODULE SOFUT - RUN EQUALS DRY OR STEP AND ',  &
     'SUBSTRUCTURE ',2A4, /33X, 'OR ONE OF THE NEW NAMES ALREADY EXISTS.')
 
 1040 FORMAT (a25,' 6220, MODULE SOFUT - RUN = GO AND SUBSTRUCTURE ',  &
     2A4,' OR ONE OF THE NEW NAMES DOES NOT EXIST')
 
 1100 RETURN
END SUBROUTINE sofut
