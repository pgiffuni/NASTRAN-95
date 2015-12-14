SUBROUTINE gp3
     
!     GP3 IS THE MAIN CONTROL PROGRAM FOR MODULE GP3.
!     IF PLOAD2 CARDS ARE PRESENT, GP3C IS EXECUTED TO BUILD PLOAD DATA
!     ON SCRATCH FILE 2 (SCR2). GP3A IS EXECUTED TO BUILD THE STATIC
!     LOADS TABLE (SLT). GP3B IS EXECUTED TO BUILD THE GRID POINT
!     TEMPERATURE TABLE (GPTT).
!     GP3D IS EXECUTED TO BUILD THE ELEMENT TEMPERATURE TABLE (ETT) FROM
!     THE GPTT AND ANY TEMPP1,TEMPP2,TEMPP3, AND TEMPRB DATA PRESENT.
 
 INTEGER :: buf1   ,buf2  ,buf  ,sysbuf,pload2,two   ,slt   ,  &
     gptt  ,geom3  ,buf3 ,STATUS,sperlk
 COMMON /BLANK / nograv ,noload,notemp
 COMMON /gp3com/ geom3 ,eqexin,geom2 ,slt   ,gptt  ,scr1  ,scr2  ,  &
     buf1  ,buf2  ,buf(50)      ,cardid(60)   ,idno(30)  &
     ,               carddt(60)   ,mask(60)     ,STATUS(60)   ,ntypes,  &
     ipload,igrav ,pload2(2)    ,load(2)      ,nopld2,  &
     temp(2)      ,tempd(2)     ,tempp1(2)    ,  &
     tempp2(2)    ,tempp3(2)    ,temprb(2)    ,buf3  , pload3(2)    ,ipld3
 COMMON /system/ sysbuf,sy(93),sperlk
 COMMON /zzzzzz/ z(1)
 COMMON /two   / two(32)
 
!     TURN PARAMETERS ON. INITIALIZE BUFFER POINTERS.
!     READ TRAILER ON GEOM3. IF PURGED, EXIT.
 
 CALL delset
 
 IF (sperlk == 0) GO TO 20
 DO  i = 1,60,2
   STATUS(i  ) =-1
   STATUS(i+1) = 0
 END DO
 20 noload = -1
 nograv = -1
 notemp = -1
 buf1   = korsz(z) - sysbuf - 2
 buf2   = buf1 - sysbuf
 buf3   = buf2 - sysbuf - 2
 buf(1) = geom3
 CALL rdtrl (buf)
 IF (buf(1) /= geom3) RETURN
 
!     IF THE SLT IS PURGED, BYPASS THE SLT PHASE OF GP3.
!     OTHERWISE, IF PLOAD2 CARDS PRESENT, EXECUTE GP3C.
!     EXECUTE GP3A TO COMPLETE SLT PHASE.
 
 buf(7) = slt
 CALL rdtrl (buf(7))
 IF (buf(7) /= slt) GO TO 30
 CALL gp3c
 CALL gp3a
 
!     IF THE GPTT IS NOT PURGED, EXECUTE GP3B TO BUILD IT.
 
 30 buf(7) = gptt
 CALL rdtrl (buf(7))
 IF (buf(7) /= gptt) RETURN
 
!     GP3B WILL FORM A GPTT ON SCR1 AND THEN GP3D WILL READ SCR1 AND
!     THE TEMPP1,TEMPP2,TEMPP3, AND TEMPRB DATA FROM GEOM3 TO FORM THE
!     ETT (ELEMENT TEMPERATURE TABLE) ON THE OUTPUT FILE GPTT.
 
 CALL gp3b
 CALL gp3d
 RETURN
END SUBROUTINE gp3
