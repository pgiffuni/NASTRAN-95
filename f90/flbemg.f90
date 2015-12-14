SUBROUTINE flbemg
     
!     GENERATES ELEMENT AREA FACTOR AND GRAVITIATIONAL STIFFNESS
!     MATRICES
 
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,mpt  &
     ,geom3    ,cstm     ,uset     ,eqexin   ,usetf  &
     ,usets    ,af       ,dkgg     ,fbelm    ,frelm  &
     ,conect   ,afmat    ,afdict   ,kgmat    ,kgdict  &
     ,z        ,grav(2)  ,pos      ,frrec(7) ,fbrec(12)  &
     ,FILE     ,NAME(2)  ,dict(2)
 
 LOGICAL :: error    ,nocard
 
 DOUBLE PRECISION :: afe(48)  ,kge(144)
 
!     GINO FILES
 
 COMMON / flbfil /       geom2    ,ect      ,bgpdt    ,sil  &
     ,mpt      ,geom3    ,cstm     ,uset ,eqexin   ,usetf    ,usets    ,af  &
     ,dkgg     ,fbelm    ,frelm    ,conect ,afmat    ,afdict   ,kgmat    ,kgdict
 
!     OPEN CORE
 
 COMMON / zzzzzz /       z(1)
 
!     CORE POINTERS
 
 COMMON / flbptr /       error    ,icore    ,lcore    ,ibgpdt  &
     ,nbgpdt   ,isil     ,nsil     ,igrav ,ngrav    ,igrid    ,ngrid    ,ibuf1  &
     ,ibuf2    ,ibuf3    ,ibuf4    ,ibuf5
 
!     MODULE PARAMETERS
 
 COMMON /BLANK/     nograv   ,nofree   ,tilt(2)
 
 DATA NAME / 4HFLBE,4HMG   /
 DATA grav / 4401 , 44 /
 
!***********************************************************************
 
!     READ MATERIAL PROPERTY DATA INTO CORE
 
 imat = icore
 nz = ibuf5 - imat
 CALL premat(z(imat),z(imat),z(ibuf1),nz,nmat,mpt,0)
 
!     READ CSTM DATA INTO CORE
 
 icstm = imat + nmat
 ncstm = 0
 nz = ibuf5 - icstm
 FILE = cstm
 CALL OPEN(*20,cstm,z(ibuf1),0)
 CALL fwdrec(*1002,cstm)
 CALL READ(*1002,*10,cstm,z(icstm),nz,0,ncstm)
 GO TO 1008
 10 CALL CLOSE(cstm,1)
 CALL pretrd(z(icstm),ncstm)
 
!     READ GRAV DATA INTO CORE
 
 20 igrav = icstm + ncstm
 ngrav = 0
 nz = ibuf5 - igrav
 nograv = -1
 nocard = .true.
 FILE = geom3
 CALL preloc(*40,z(ibuf1),geom3)
 CALL locate(*30,z(ibuf1),grav,id)
 nocard = .false.
 CALL READ(*1002,*30,geom3,z(igrav),nz,0,ngrav)
 GO TO 1008
 
 30 CALL CLOSE(geom3,1)
 40 CONTINUE
 
!     OPEN MATRIX AND DICTIONARY FILES
 
 CALL gopen(afmat,z(ibuf2),1)
 CALL gopen(afdict,z(ibuf4),1)
 IF(nocard) GO TO 60
 CALL gopen(kgmat,z(ibuf3),1)
 CALL gopen(kgdict,z(ibuf5),1)
 
 
!     PASS THROUGH FBELM FILE AND PROCESS EACH ENTRY ON THE BOUNDARY.
!     SUBROUTINE BOUND WILL GENERATE THE ELEMENT MATRICES FOR
!     EACH ENTRY.
 
 60 FILE = fbelm
 CALL gopen(fbelm,z(ibuf1),0)
 70 CALL READ(*1002,*120,fbelm,fbrec,12,0,n)
 
 CALL bound(fbrec,afe,nafe,kge,nkge)
 IF(error) GO TO 70
 
!     CONVERT GRID POINTS TO SILS
 
 DO  i=1,4
   j = fbrec(i+2) - 1
   IF(j >= 0) fbrec(i+2) = z(isil+j)
   j = fbrec(i+8) - 1
   IF(j >= 0) fbrec(i+8) = z(isil+j)
 END DO
 
!     WRITE AREA MATRICES AND DICTIONARY ENTRUES
 
 CALL WRITE(afmat,fbrec(3),4,0)
 CALL WRITE(afmat,fbrec(9),4,0)
 CALL WRITE(afmat,afe,nafe,1)
 CALL savpos(afmat,pos)
 dict(2) = pos
 DO  i=1,4
   dict(1) = fbrec(i+8)
   IF(dict(1) < 0) CYCLE
   CALL WRITE(afdict,dict,2,0)
 END DO
 
!     WRITE GRAVITATIONAL STIFFNESS MATRICES IF THEY EXIST
 
 IF(nkge == 0) GO TO 70
 CALL WRITE(kgmat,fbrec(3),4,0)
 CALL WRITE(kgmat,fbrec(3),4,0)
 CALL WRITE(kgmat,kge,nkge,1)
 CALL savpos(kgmat,pos)
 dict(2) = pos
 DO  i=1,4
   jsil = fbrec(i+2)
   IF(jsil < 0) CYCLE
   DO  j=1,3
     dict(1) = jsil
     CALL WRITE(kgdict,dict,2,0)
     jsil = jsil + 1
   END DO
 END DO
 
 GO TO 70
 120 CALL CLOSE(fbelm,1)
 
 
!     PASS THROUGH FRELM FILE AND PROCESS EACH ENTRY ON THE FREE
!     SURFACE.  SUBROUTINE FLFREE WILL CALCULATE THE AREA AND
!     GRAVITATIONAL STIFFNESS MATRICES FOR EACH ENTRY
 
 IF(nofree < 0) GO TO 180
 FILE = frelm
 CALL gopen(frelm,z(ibuf1),0)
 130 CALL READ(*1002,*170,frelm,frrec,7,0,n)
 
 CALL flfree(frrec,afe,nafe,kge,nkge)
 IF(error) GO TO 130
 
!     CONVERT GRID POINTS TO SILS
 
 DO  i=1,4
   j = frrec(i+2) - 1
   IF(j >= 0) frrec(i+2) = z(isil+j)
 END DO
 
!     WRITE AREA MATRICES AND DICTIONARY ENTRIES
 
 CALL WRITE(afmat,frrec(3),4,0)
 CALL WRITE(afmat,frrec(3),4,0)
 CALL WRITE(afmat,afe,nafe,1)
 CALL savpos(afmat,pos)
 dict(2) = pos
 DO  i=1,4
   dict(1) = frrec(i+2)
   IF(dict(1) < 0) CYCLE
   CALL WRITE(afdict,dict,2,0)
 END DO
 
!     WRITE GRAVITATIONAL STIFFNESS MATRICES IF THEY EXIST
 
 IF(nkge == 0) GO TO 130
 CALL WRITE(kgmat,frrec(3),4,0)
 CALL WRITE(kgmat,frrec(3),4,0)
 CALL WRITE(kgmat,kge,nkge,1)
 CALL savpos(kgmat,pos)
 dict(2) = pos
 DO  i=1,4
   dict(1) = frrec(i+2)
   IF(dict(1) < 0) CYCLE
   CALL WRITE(kgdict,dict,2,0)
 END DO
 
 GO TO 130
 170 CALL CLOSE(frelm,1)
 
!     CLOSE FILES AND RETURN
 
 180 CALL CLOSE(afmat,1)
 CALL CLOSE(afdict,1)
 IF(nocard) GO TO 190
 CALL CLOSE(kgmat,1)
 CALL CLOSE(kgdict,1)
 
 190 CONTINUE
 RETURN
 
!     ERROR CONDITIONS
 
 1002 n = -2
 GO TO 1100
 1008 n = -8
 1100 CALL mesage(n,FILE,NAME)
 RETURN
END SUBROUTINE flbemg
