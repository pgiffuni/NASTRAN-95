SUBROUTINE dpd
     
!     DPD IS MAIN CONTROL PROGRAM FOR THE DYNAMICS POOL DISTRIBUTOR.
 
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool ,  &
     dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  ,  &
     scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  ,  &
     buf4  ,sysbuf,ngrid ,epoint,seqep ,z     ,loads ,  &
     eqdyn ,dload ,freq1 ,freq  ,tic   ,tstep ,tf    ,  &
     psd   ,eigr  ,eigb  ,eigc  ,sdt
 DIMENSION       buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)   ,  &
     nam(2)    ,loads(32)    ,dload(2)     ,freq1(2) ,  &
     freq(2)   ,zz(1)        ,bufr(20)     ,nolin(21),  &
     tic(2)    ,tstep(2)     ,tf(2)        ,psd(2)   ,  &
     msg(3)    ,eigr(2)      ,eigb(2)      ,eigc(2)
 COMMON /BLANK / luset ,lusetd,notfl ,nodlt ,nopsdl,nofrl ,nonlft,  &
     notrl ,noeed ,nosdt ,noue
 COMMON /dpdcom/ dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd ,  &
     dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed   ,  &
     scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,  &
     buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn,  &
     loads ,dload ,freq1 ,freq  ,nolin ,nogo  ,  &
     msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  ,  &
     eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineqc
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf
 EQUIVALENCE     (z(1),zz(1)),(buf(1),bufr(1)),(msg(2),ngrid)
 
!     INITIALIZE CONTROL PARAMETERS.
 
 notfl  = -1
 nodlt  = -1
 nopsdl = -1
 nofrl  = -1
 nonlft = -1
 notrl  = -1
 noeed  = -1
 nosdt  = -1
 noue   = -1
 nogo   =  0
 ineq   =  0
 DO  i = 1,7
   mcb(i) = 0
 END DO
 
!     PERFORM BUFFER ALLOCATION
 
 buf1 = korsz(z) - sysbuf - 2
 buf2 = buf1 - sysbuf
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 
!     IF DYNAMICS POOL IS PURGED, EXIT. OTHERWISE, EXECUTE THE PHASES
!     OF DPD
 
 buf(1) = dpool
 CALL rdtrl (buf)
 IF (buf(1) /= dpool) RETURN
 CALL dpd1
 CALL dpd2
 CALL dpd3
 CALL dpd4
 CALL dpd5
 IF (nogo /= 0) CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE dpd
