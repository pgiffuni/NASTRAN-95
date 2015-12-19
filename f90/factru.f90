SUBROUTINE factru(*,a,lll,ull,scr1,scr2,scr3)
     
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN)                      :: lll
 INTEGER, INTENT(IN)                      :: ull
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 
 DOUBLE PRECISION :: dett,mindia
 
 COMMON /dcompx/ia(7),il(7),iu(7),iscr1,iscr2,iscr3,dett,ipow,  &
     nz,mindia,ib,ibb
 COMMON /system/sys(54),iprec
 COMMON /zzzzzz/ xx(1)
 
! ----------------------------------------------------------------------
 
 ib = 0
 ibb = 0
 ia(1) = a
 CALL rdtrl(ia)
 il(1)=lll
 iu(1)=ull
 iscr1=scr1
 iscr2=scr2
 iscr3=scr3
 nz = korsz(xx)
 il(3) = ia(3)
 iu(3) = ia(3)
 il(4) =4
 iu(4) =5
 iu(5) = iprec
 il(5) = iprec
 CALL decomp(*10,xx,xx,xx)
 CALL wrttrl(il)
 CALL wrttrl(iu)
 RETURN
 10 RETURN 1
 
END SUBROUTINE factru
