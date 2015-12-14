SUBROUTINE copy
     
!     COPY  INPUT /OUTPUT/ PARAM $
 
!     THIS UTILITY MODULE GENERATES A PHYSICAL COPY OF THE INPUT DATA
!     BLOCK IF THE VALUE OF PARAM IS LESS THAN ZERO (DEFAULT IS -1).
!     THE OUTPUT DATA BLOCK CARRIES THE INPUT DATA BLOCK NAME IN THE
!     HEADER RECORD.
!     IF PARAM IS SET TO ZERO, THE OUTPUT DATA BLOCK WILL HAVE ITS OWN
!     NAME IN THE OUTPUT FILE HEADER RECORD.  (IMPLEMENTED IN JUNE 84)
 
 
 INTEGER :: modnam(2),sysbuf,output,itrl(7),in(15),out(15)
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / iparam
 COMMON /xfist / ifist(1)
 COMMON /xfiat / ifiat(1)
 DATA    INPUT / 101 /, output / 201 /, modnam / 4HCOPY,4H    /
 
!     RETURN IF IPARAM NOT GREATER THAN ZERO
 
 IF (iparam == 0) iparam = -1111
 IF (iparam >= 0) RETURN
 
!     COMPUTE OPEN CORE AND INITIALIZE GINO BUFFERS
 
 nzwd = korsz(z(1))
 IF (nzwd <= 0) CALL mesage (-8,0,modnam)
 ibuf1 = nzwd  - sysbuf
 ibuf2 = ibuf1 - sysbuf
 lcore = ibuf2 - 1
 IF (lcore <= 0) CALL mesage (-8,0,modnam)
 
!     OPEN INPUT AND OUTPUT DATA BLOCKS
 
 in(1)   = INPUT
 out(1)  = output
 itrl(1) = 101
 CALL rdtrl (itrl)
 CALL OPEN  (*1001,INPUT,z(ibuf1),0)
 CALL OPEN  (*1002,output,z(ibuf2),1)
 CALL cpyfil (in,out,z(1),lcore,icount)
 CALL CLOSE (output,1)
 CALL CLOSE (INPUT,1)
 itrl(1) = 201
 CALL wrttrl (itrl)
 RETURN
 
 1001 CALL mesage (-1,INPUT,modnam)
 1002 CALL mesage (-1,output,modnam)
 RETURN
END SUBROUTINE copy
