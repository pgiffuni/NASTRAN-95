SUBROUTINE gpwg
     
!     GRID POINT WEIGHT GENERATOR
 
!     INPUTS  - BGPDT,CSTM,EQEXIN,MGG
 
!     OUTPUTS - OGPWG
 
!     PARAMETERS -- POINT,WTMASS
 
 INTEGER :: bgpdt,cstm,eqexin,ogpwg,scr1,scr2,scr3,scr4,point
 COMMON /BLANK/ point,wtmass
 DATA    bgpdt, cstm,eqexin,mgg, ogpwg, scr1,scr2,scr3,scr4 /  &
     101  , 102 ,103   ,104, 201  , 301 ,302 ,303 ,304  /
 
!     FORM D MATRIX (TRANSPOSED)
 
 ip = point
 
 CALL gpwg1a (point,bgpdt,cstm,eqexin,scr3,nogo)
 
!     CHECK FOR AN ALL SCALAR PROBLEM AND A STUPID USER
 
 IF (nogo == 0) GO TO 10
 
!     COMPUTE MZERO = DT*MGG*D
 
 CALL tranp1 (scr3,scr1,2,scr2,scr4,0,0,0,0,0,0)
 CALL ssg2b  (mgg ,scr1,0,scr2,0,1,1,scr3)
 CALL ssg2b  (scr1,scr2,0,scr4,1,1,1,scr3)
 
!     M-ZERO IS ON SCR4
 
!     FORM OUTPUT  STUFF
 
 IF (point == 0) ip = 0
 CALL gpwg1b (scr4,ogpwg,wtmass,ip)
 10 RETURN
END SUBROUTINE gpwg
