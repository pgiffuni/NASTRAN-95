SUBROUTINE ssgetd (elid,ti,grids)
     
!     THIS ROUTINE (CALLED BY -EDTL-) READS ELEMENT TEMPERATURE
!     DATA FROM A PRE-POSITIONED RECORD
 
!     ELID   = ID OF ELEMENT FOR WHICH DATA IS DESIRED
!     TI     = BUFFER DATA IS TO BE RETURNED IN
!     GRIDS  = 0 IF EL-TEMP FORMAT DATA IS TO BE RETURNED
!            = NO. OF GRID POINTS IF GRID POINT DATA IS TO BE RETURNED.
!     ELTYPE = ELEMENT TYPE TO WHICH -ELID- BELONGS
!     OLDEL  = ELEMENT TYPE CURRENTLY BEING WORKED ON (INITIALLY 0)
!     EORFLG = .TRUE. WHEN ALL DATA HAS BEEN EXHAUSTED IN RECORD
!     ENDID  = .TRUE. WHEN ALL DATA HAS BEEN EXHAUSTED WITHIN AN ELEMENT
!              TYPE.
!     BUFFLG = NOT USED
!     ITEMP  = TEMPERATURE LOAD SET ID
!     IDEFT  = NOT USED
!     IDEFM  = NOT USED
!     RECORD = .TRUE. IF A RECORD OF DATA IS INITIALLY AVAILABLE
!     DEFALT = THE DEFALT TEMPERATURE VALUE OR -1 IF IT DOES NOT EXIST
!     AVRAGE = THE AVERAGE ELEMENT TEMPERATURE
 
 
 INTEGER, INTENT(IN OUT)                  :: elid
 INTEGER, INTENT(OUT)                     :: ti(7)
 INTEGER, INTENT(IN)                      :: grids
 LOGICAL :: eorflg   ,endid    ,bufflg   ,record
 INTEGER :: eltype   ,oldel   , NAME(2)  ,gptt     ,defalt
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ dum      ,iout
 COMMON /ssgett/ eltype   ,oldel    ,eorflg   ,endid    ,bufflg  ,  &
     itemp    ,ideft    ,idefm    ,record
 COMMON /loadx / dummy(9) ,gptt
 COMMON /fpt   / defalt
 DATA    NAME  / 4HSSGE,4HTD   /    ,maxwds   / 15 /
 
 IF (itemp /= 0) GO TO 20
 DO  i = 1,maxwds
   ti(i) = 0
 END DO
 RETURN
 
 20 IF (.NOT.record .OR. eorflg) GO TO 80
 30 IF (eltype /= oldel) GO TO 150
 IF (endid) GO TO 80
 
!     HERE WHEN ELTYPE IS AT HAND AND END OF THIS TYPE DATA
!     HAS NOT YET BEEN REACHED.  READ AN ELEMENT ID
 
 40 CALL READ (*200,*210,gptt,id,1,0,flag)
 IF (id == 0) THEN
   GO TO    80
 END IF
 50 IF (IABS(id) == elid) IF (id) 90,90,70
 IF (id > 0) THEN
   GO TO    60
 ELSE
   GO TO    40
 END IF
 60 CALL READ (*200,*210,gptt,ti,nwords,0,flag)
 GO TO 40
 
!     MATCH ON ELEMNT ID MADE, AND IT WAS WITH DATA.
!     IF QUAD4 OR TRIA3 ELEMENT, SET THE TI(7) FLAG FOR TLQD4D/S (QAUD4)
!     OR TLTR3D/S (TRIA3)
 
 70 CALL READ (*200,*210,gptt,ti,nwords,0,flag)
 IF (eltype /= 64 .AND. eltype /= 83) RETURN
 ti(7) = 13
 IF (ti(6) /= 1) ti(7) = 2
 RETURN
 
!     NO MORE DATA FOR THIS ELEMENT TYPE
 
 80 endid = .true.
 
!     NO DATA FOR ELEMENT ID DESIRED, THUS USE DEFALT
 
 90 IF (defalt == -1) GO TO 130
 IF (grids  >  0) GO TO 110
 DO  i = 2,maxwds
   ti(i) = 0
 END DO
 ti(1) = defalt
 IF (eltype == 34) ti(2) = defalt
 RETURN
 110 DO  i = 1,grids
   ti(i) = defalt
 END DO
 ti(grids+1) = defalt
 RETURN
 
!     NO TEMP DATA OR DEFALT
 
 130 WRITE  (iout,140) ufm,elid,itemp
 140 FORMAT (a23,' 4017. THERE IS NO TEMPERATURE DATA FOR ELEMENT',i9,  &
     ' IN SET',i9)
 CALL mesage (-61,0,0)
 
!     LOOK FOR MATCH ON ELTYPE (FIRST SKIP ANY UNUSED ELEMENT DATA)
 
 150 IF (endid) GO TO 180
 160 CALL READ (*200,*210,gptt,id,1,0,flag)
 IF (id < 0) THEN
   GO TO   160
 ELSE IF (id == 0) THEN
   GO TO   180
 END IF
 170 CALL READ (*200,*210,gptt,ti,nwords,0,flag)
 GO TO 160
 
!     READ ELTYPE AND COUNT
 
 180 CALL READ (*200,*190,gptt,ti,2,0,flag)
 oldel  = ti(1)
 nwords = ti(2)
 endid  = .false.
 GO TO 30
 
!     END OF RECORD HIT
 
 190 eorflg = .true.
 GO TO 80
 
 200 CALL mesage (-2,gptt,NAME)
 210 CALL mesage (-3,gptt,NAME)
 RETURN
END SUBROUTINE ssgetd
