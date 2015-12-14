SUBROUTINE ds1etd (elid,ti,grids)
     
!     THIS ROUTINE (CALLED BY -DS1-) READS ELEMENT TEMPERATURE
!     DATA FROM A PRE-POSITIONED RECORD
 
!     ELID   = ID OF ELEMENT FOR WHICH DATA IS DESIRED
!     TI     = BUFFER DATA IS TO BE RETURNED IN
!     GRIDS  = 0 IF EL-TEMP FORMAT DATA IS TO BE RETURNED
!            = NO. OF GRID POINTS IF GRID POINT DATA IS TO BE RETURNED.
!     ELTYPE = ELEMENT TYPE TO WHICH -ELID- BELONGS
!     OLDEL  = ELEMENT TYPE CURRENTLY BEING WORKED ON (INITIALLY 0)
!     OLDEID = ELEMENT ID FROM LAST CALL
!     EORFLG =.TRUE. WHEN ALL DATA HAS BEEN EXHAUSTED IN RECORD
!     ENDID  =.TRUE. WHEN ALL DATA HAS BEEN EXHAUSTED WITHIN AN ELEMENT
!              TYPE.
!     BUFFLG = NOT USED
!     ITEMP  = TEMPERATURE LOAD SET ID
!     IDEFT  = NOT USED
!     IDEFM  = NOT USED
!     RECORD =.TRUE. IF A RECORD OF DATA IS INITIALLY AVAILABLE
!     DEFALT = THE DEFALT TEMPERATURE VALUE OR -1 IF IT DOES NOT EXIST
!     AVRAGE = THE AVERAGE ELEMENT TEMPERATURE
 
 
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(OUT)                     :: ti(2)
 INTEGER, INTENT(IN)                      :: grids
 LOGICAL :: eorflg   ,endid    ,bufflg   ,record
 INTEGER :: oldeid   , eltype   , oldel    ,NAME(2)  ,gptt     ,defalt
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ dum      ,iout
 COMMON /ds1ett/ eltype   ,oldel    ,eorflg   ,endid    ,bufflg   ,  &
     itemp    ,defalt   ,iback    ,record   ,oldeid
 DATA    NAME  / 4HDS1E,4HTD  /,  maxwds / 33 /,  gptt  / 102 /
 
 IF (oldeid == elid) RETURN
 oldeid = elid
 
 IF (itemp > 0) GO TO 20
 DO  i = 1,maxwds
   ti(i) =-1
 END DO
 RETURN
 
 20 IF (.NOT.record .OR. eorflg) GO TO 50
 15 IF (eltype /= oldel) GO TO 30
 IF (endid) GO TO 50
 
!     HERE WHEN ELTYPE IS AT HAND AND END OF THIS TYPE DATA
!     HAS NOT YET BEEN REACHED.  READ AN ELEMENT ID
 
 35 CALL READ (*5002,*5001,gptt,id,1,0,flag)
 IF (id == 0) THEN
   GO TO    50
 END IF
 40 IF (IABS(id) == elid) IF (id) 51,51,70
 IF (id > 0) THEN
   GO TO    45
 ELSE
   GO TO    35
 END IF
 45 CALL READ (*5002,*5001,gptt,ti,nwords,0,flag)
 GO TO 35
 
!     MATCH ON ELEMNT ID MADE AND IT WAS WITH DATA
 
 70 CALL READ (*5002,*5001,gptt,ti,nwords,0,flag)
 RETURN
 
!     NO MORE DATA FOR THIS ELEMENT TYPE
 
 50 endid = .true.
 
!     NO DATA FOR ELEMENT ID DESIRED, THUS USE DEFALT
 
 51 IF (defalt == -1) GO TO 100
 IF (grids  >  0) GO TO 75
 DO  i = 2,maxwds
   ti(i) = 0
 END DO
 ti(1) = defalt
 IF (eltype == 34) ti(2) = defalt
 RETURN
 
 75 DO  i = 1,grids
   ti(i) = defalt
 END DO
 ti(grids+1) = defalt
 RETURN
 
!     NO TEMP DATA OR DEFALT
 
 100 WRITE  (iout,301) ufm,elid,itemp
 301 FORMAT (a23,' 4016, THERE IS NO TEMPERATURE DATA FOR ELEMENT',i9,  &
     ' IN SET',i9)
 CALL mesage (-61,0,0)
 
!     LOOK FOR MATCH ON ELTYPE (FIRST SKIP ANY UNUSED ELEMENT DATA)
 
 30 IF (endid) GO TO 32
 31 CALL READ (*5002,*5001,gptt,id,1,0,flag)
 IF (id < 0) THEN
   GO TO    31
 ELSE IF (id == 0) THEN
   GO TO    32
 END IF
 33 CALL READ (*5002,*5001,gptt,ti,nwords,0,flag)
 GO TO 31
 
!     READ ELTYPE AND COUNT
 
 32 CALL READ (*5002,*300,gptt,ti,2,0,flag)
 oldel  = ti(1)
 nwords = ti(2)
 endid  = .false.
 iback  = 1
 GO TO 15
 
!     END OF RECORD HIT
 
 300 eorflg = .true.
 GO TO 50
 5002 CALL mesage (-2,gptt,NAME)
 5001 CALL mesage (-3,gptt,NAME)
 RETURN
END SUBROUTINE ds1etd
