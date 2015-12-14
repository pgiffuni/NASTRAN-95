SUBROUTINE sdretd (elid,ti,grids)
     
!     THIS ROUTINE (CALLED BY -SDR2E-) READS ELEMENT TEMPERATURE
!     DATA FROM A PRE-POSITIONED RECORD
 
!     ELID   = ID OF ELEMENT FOR WHICH DATA IS DESIRED
!     TI     = BUFFER DATA IS TO BE RETURNED IN
!     GRIDS  = 0 IF EL-TEMP FORMAT DATA IS TO BE RETURNED
!            = NO. OF GRID POINTS IF GRID POINT DATA IS TO BE RETURNED.
!     ELTYPE = ELEMENT TYPE TO WHICH -ELID- BELONGS
!     OLDEL  = ELEMENT TYPE CURRENTLY BEING WORKED ON (INITIALLY 0)
!     OLDEID = ELEMENT ID FROM LAST CALL
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
 
 
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(OUT)                     :: ti(33)
 INTEGER, INTENT(IN)                      :: grids
 LOGICAL :: eorflg   ,endid    ,bufflg   ,record
 INTEGER :: oldeid   , eltype   , oldel    ,NAME(2)  ,gptt     ,defalt
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ dum      ,iout
 COMMON /sdrett/ eltype   ,oldel    ,eorflg   ,endid    ,bufflg   ,  &
     itemp    ,ideft    ,idefm    ,record   ,oldeid
 COMMON /sdr2x2/ dummy(6) ,gptt     ,dum20(20)
 COMMON /sdr2de/ dum2(2)  ,defalt
 DATA    NAME  / 4HSDRE   ,4HTD  /  ,maxwds   / 33 /
 
 IF (oldeid == elid) RETURN
 oldeid = elid
 IF (itemp /= 0) GO TO 20
 DO  i = 1,maxwds
   ti(i) = 0
 END DO
 RETURN
 
 20 IF (.NOT.record .OR. eorflg) GO TO 80
 IF (eltype /= oldel) GO TO 160
 IF (endid) GO TO 80
 
!     HERE WHEN ELTYPE IS AT HAND AND END OF THIS TYPE DATA
!     HAS NOT YET BEEN REACHED.  READ AN ELEMENT ID
 
 40 CALL READ (*300,*310,gptt,id,1,0,flag)
 IF (id == 0) THEN
   GO TO    80
 END IF
 50 IF (IABS(id) == elid) IF (id) 90,90,70
 IF (id > 0) THEN
   GO TO    60
 ELSE
   GO TO    40
 END IF
 60 CALL READ (*300,*310,gptt,ti,nwords,0,flag)
 GO TO 40
 
!     MATCH ON ELEMNT ID MADE AND IT WAS WITH DATA
 
 70 CALL READ (*300,*310,gptt,ti,nwords,0,flag)
 
!     IF QUAD4 (ELTYPE 64) OR TRIA3 (ELTYPE 83) ELEMENT, SET FLAG FOR
!     SQUD42 OR STRI32
 
 IF (eltype /= 64 .OR. eltype /= 83) RETURN
 ti(7) = 13
 IF (ti(6) /= 1) ti(7) = 2
 RETURN
 
!     NO MORE DATA FOR THIS ELEMENT TYPE
 
 80 endid = .true.
 
!     NO DATA FOR ELEMENT ID DESIRED, THUS USE DEFALT
 
 90 IF (defalt == -1) GO TO 140
 IF (grids  >  0) GO TO 110
 DO  i = 2,maxwds
   ti(i) = 0
 END DO
 ti(1) = defalt
 IF (eltype == 34) ti(2) = defalt
 RETURN
 
 110 IF (eltype /= 64 .OR. eltype /= 83) GO TO 120
!                 QUAD4             TRIA3
 ti(4) = 0
 ti(5) = 0
 ti(6) = 0
 ti(7) = 0
 120 DO  i = 1,grids
   ti(i) = defalt
 END DO
 ti(grids+1) = defalt
 RETURN
 
!     NO TEMP DATA OR DEFALT
 
 140 WRITE  (iout,150) ufm,elid,itemp
 150 FORMAT (a23,' 4016, THERE IS NO TEMPERATURE DATA FOR ELEMENT',i9,  &
     ' IN SET',i9)
 CALL mesage (-61,0,0)
 
!     LOOK FOR MATCH ON ELTYPE (FIRST SKIP ANY UNUSED ELEMENT DATA)
 
 160 IF (endid) GO TO 190
 170 CALL READ (*300,*310,gptt,id,1,0,flag)
 IF (id < 0) THEN
   GO TO   170
 ELSE IF (id == 0) THEN
   GO TO   190
 END IF
 180 CALL READ (*300,*310,gptt,ti,nwords,0,flag)
 GO TO 170
 
!     READ ELTYPE AND COUNT
 
 190 CALL READ (*300,*200,gptt,ti,2,0,flag)
 oldel  = ti(1)
 nwords = ti(2)
 endid = .false.
 GO TO 40
!     END OF RECORD HIT
 
 200 eorflg = .true.
 GO TO 80
 
 300 CALL mesage (-2,gptt,NAME)
 310 CALL mesage (-3,gptt,NAME)
 RETURN
END SUBROUTINE sdretd
