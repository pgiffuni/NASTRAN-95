SUBROUTINE sdhtf2(ieqex,neqex)
!*****
!     THIS ROUTINE CALCULATES TEMPERATURE GRADIENTS AND HEAT FLOWS
!     FOR ALL ELEMENTS IN A HEAT TRANSFER PROBLEM.
!      DATA IS OUTPUT FOR ELEMENT FORCE REQUEST ONLY.
!******
 
 INTEGER, INTENT(IN)                      :: ieqex
 INTEGER, INTENT(IN)                      :: neqex
 INTEGER :: igrad(3), iqout(3), ftube
 REAL :: esta(202)
 DIMENSION iz(1),ipt(21)
 COMMON /zzzzzz/  zz(1)
 COMMON  /sdr2x4/  dummy(35),ivec
 COMMON/sdr2x7/ide,isil(32),nq,nsil,NAME(2),rk(9),ce(96),  &
     dum(58),ido,namo(2),tgrad(3),qout(3)
 COMMON/sdr2x8/tvec(32)
 EQUIVALENCE  (tgrad(1),igrad(1)) ,(qout(1),iqout(1))
 EQUIVALENCE  (zz(1),iz(1)), (esta(1),ide)
 DATA ihex/4HIHEX/,ione,itwo,ithr/4H1   ,4H2   ,4H3   /
 DATA ihex1,ihex2,ihex3/4HHEX1,4HHEX2,4HHEX3/
 DATA ftube/4HFTUB/
 DATA iold/0/
 DATA ipt/4H   1,4H  e1,4H   4,4H  e2,4H   7,4H  e3,4H  10,  &
     4H  e4,4H  e5,4H  e6,4H  e7,4H  e8,4H  21,4H  e9,  &
     4H  24,4H e10,4H  27,4H e11,4H  30,4H e12,4H   0/
 
 IF (NAME(1) == ftube) GO TO 70
 DO  i=1,3
   igrad(i)= 1
   iqout(i)= 1
 END DO
 ido= ide
 namo(1)= NAME(1)
 namo(2)= NAME(2)
 
! FOR ISOPARAMETRIC SOLIDS, GET SIL NUMBER AND CONVERT TO EXTERNAL.
! STORE IT IN NAMO(2)
 
 IF(namo(1) /= ihex) GO TO 29
 IF(iold == ide) GO TO 11
 iold=ide
 istrpt=0
 11 IF(namo(2) == ione) namo(1)=ihex1
 IF(namo(2) == itwo) namo(1)=ihex2
 IF(namo(2) == ithr) namo(1)=ihex3
 istrpt=istrpt+1
 IF(istrpt == nsil+1.OR.istrpt == 21) iold=0
 IF(namo(1) == ihex3) GO TO 12
 IF(namo(1) == ihex1.AND.istrpt == 9) GO TO 15
 IF(namo(1) == ihex2.AND.istrpt == 21) GO TO 15
 GO TO 13
 12 namo(2)=ipt(istrpt)
 GO TO 29
 13 isub1=ieqex+1
 isub2=ieqex+neqex-1
 DO  jjj=isub1,isub2,2
   ns=iz(jjj)/10
   IF(ns /= isil(istrpt)) CYCLE
   namo(2)=iz(jjj-1)
   GO TO 29
 END DO
 CALL mesage(-30,164,iz(jjj))
 15 namo(2)=0
 29 CONTINUE
 IF(nq <= 0) GO TO 60
 DO  i=1,nsil
   tvec(i)= 0.0
   ip= isil(i)
   IF( ip == 0) CYCLE
   itemp = ivec + ip -1
   tvec(i) = zz(itemp)
 END DO
!***
 CALL gmmats( ce(1),nq,nsil,0, tvec(1),nsil,1,0, tgrad(1) )
 
 CALL gmmats( rk(1),nq,nq,0, tgrad(1),nq,1,0, qout(1) )
 
 DO  i=1,nq
   qout(i) =-qout(i)
 END DO
 RETURN
 60 tgrad(1) = 0.0
 qout(1) = 0.0
 GO TO 80
 
 70 ido=ide
 itemp=ivec + isil(1) - 1
 tvec(1)=zz(itemp)
 esta(202)=tvec(1)*esta(4)
 
 80 RETURN
END SUBROUTINE sdhtf2
