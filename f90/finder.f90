SUBROUTINE finder( nam , subno , comno )
     
 
!     THIS SUBROUTINE READS THE TABLE OF CONTENTS OF SUBSTRUCTURES
!     BEING COMBINED ( SCRATCH FILE SCTOC ) AND FOR ANY GIVEN
!     BASIC SUBSTRUCTURE NAME ( NAM ) RETURNS THE ID NUMBER OF THE
!     PSEUDO-STRUCTURE CONTAINING IT ( SUBNO ) AND ITS POSITION IN
!     THE COMPONENT LIST FOR THAT STRUCTURE ( COMNO ).  IF A NAME
!     DOES NOT APPEAR IN THE SCTOC AN ERROR MESSAGE IS ISSUED.
 
 
 INTEGER, INTENT(IN)                      :: nam(2)
 INTEGER, INTENT(OUT)                     :: subno
 INTEGER, INTENT(OUT)                     :: comno
 INTEGER :: sctoc,buf4,id(3), cnam(2),outt
 LOGICAL :: tocopn
 COMMON/cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon, sctoc,geom4,casecc
 COMMON/zzzzzz/ z(1)
 COMMON/cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON/cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint,tocopn
 COMMON/cmbfnd/ inam(2),ierr
 
!     OPEN SCTOC FILE
 
 ierr = 0
 IF(.NOT.tocopn)CALL OPEN(*2001,sctoc,z(buf4),0)
 CALL REWIND( sctoc )
 
 loop1:  DO  i=1,npsub
   CALL READ(*2001,*2002,sctoc,id,3,0,nnn)
   ncom = id(3)
   DO  j=1,ncom
     IEOR = 0
     IF( j == ncom ) IEOR = 1
     CALL READ(*2001,*2002,sctoc,cnam,2,IEOR,nnn)
     IF( nam(1) == cnam(1) .AND. nam(2) == cnam(2) ) GO TO 11
   END DO
 END DO loop1
 
!     IERR = 1 MEANS THAT THE SUBSTRUCTURE NAME IS NOT IN THE TOC
 
 ierr = 1
 RETURN
 11    subno = i
 inam(1) = id(1)
 inam(2) = id(2)
 comno = j
 IF( .NOT. tocopn ) CALL CLOSE( sctoc , 1 )
 2001  CONTINUE
 2002  CONTINUE
 RETURN
END SUBROUTINE finder
