SUBROUTINE bdat05
     
!     THIS SUBROUTINE PROCESSES THE GNEW BULK DATA
 
 LOGICAL :: tdat
 INTEGER :: scr2,buf3,scbdat,buf2,buf1,gnew(2),flag,geom4,  &
     conset,aaa(2),z,outt,score
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ z(1)
 COMMON /cmb001/ scr1,scr2,scbdat,scsfil,scconn,scmcon, sctoc,geom4,casecc
 COMMON /cmb002/ buf1,buf2,buf3,buf4,buf5,score,lcore,inpt,outt
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub
 COMMON /cmb004/ tdat(6)
 COMMON /BLANK / step,idry
 DATA    gnew  / 1410,14 / , aaa/ 4HBDAT,4H05   /
 
 ifile = scr2
 CALL OPEN (*30,scr2,z(buf3),1)
 ifile = geom4
 CALL locate (*20,z(buf1),gnew,flag)
 WRITE  (outt,10) ufm
 10   FORMAT (a23,' 6532, THE GNEW OPTION IS NOT CURRENTLY AVAILABLE.')
 idry = -2
 RETURN
 
 20   CALL eof (scbdat)
 CALL CLOSE (scr2,1)
 RETURN
 
 30   imsg = -1
 CALL mesage (imsg,ifile,aaa)
 RETURN
END SUBROUTINE bdat05
