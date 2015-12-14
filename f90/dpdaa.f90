SUBROUTINE dpdaa
!*****
! DPDAA PERFORMS A BINARY SEARCH IN EQDYN AND CONVERTS THE GRID NO
! AND COMPONENT CODE TO AN SIL VALUE.
!*****
 
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool  &
     ,dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  &
     ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  &
     ,buf4  ,ngrid ,eqdyn ,epoint,seqep ,z     ,loads  &
     ,psd   ,dload ,freq1 ,freq  ,tic   ,tstep ,tf ,eigr  ,eigb  ,eigc
 
 DIMENSION buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)  &
     ,nam(2)    ,loads(32)    ,dload(2)     ,freq1(2)  &
     ,freq(2)   ,zz(1)        ,bufr(20)     ,nolin(21)  &
     ,tic(2)    ,tstep(2)     ,tf(2)        ,psd(2)  &
     ,msg(3)    ,eigr(2)      ,eigb(2)      ,eigc(2)
 
 COMMON/dpdcom/dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd  &
     ,dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed  &
     ,scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  &
     ,buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn  &
     ,loads ,dload ,freq1 ,freq  ,nolin ,nogo  &
     ,msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  &
     ,eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineq
 
 COMMON/zzzzzz/z(1)
 
 EQUIVALENCE   (z(1) ,zz(1)),(buf(1),bufr(1)),(msg(2),ngrid)
 
!*****
! IF EQDYN IS NOT IN CORE, READ IT IN AND SET FLAG.
!*****
 IF(ineq /= 0) GO TO 1
 CALL gopen(eqdyn,z(buf3),0)
 CALL fread(eqdyn,z,neqdyn+1,1)
 CALL CLOSE(eqdyn,1)
 ineq= 1
!*****
! PERFORM SEARCH.
!*****
 1 klo= 1
 khi= kn
 ngrid= buf(l)
 2 k= (klo+khi+1)/2
 3 IF(ngrid - z(2*k-1) < 0) THEN
   GO TO     4
 ELSE IF (ngrid - z(2*k-1) == 0) THEN
   GO TO    11
 ELSE
   GO TO     5
 END IF
 4 khi= k
 GO TO 6
 5 klo= k
 6 IF(khi-klo-1 < 0) THEN
   GO TO    10
 ELSE IF (khi-klo-1 == 0) THEN
   GO TO     7
 ELSE
   GO TO     2
 END IF
 7 IF(k == klo) GO TO 8
 k= klo
 GO TO 9
 8 k= khi
 9 klo= khi
 GO TO 3
 10 CALL mesage(30,msg,msg(2))
 nogo= 1
 11 buf(l)= z(2*k)
 IF(buf(l+1) /= 0) buf(l)= buf(l)+buf(l+1)-1
 RETURN
END SUBROUTINE dpdaa
