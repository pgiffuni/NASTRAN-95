SUBROUTINE rand1(FILE,mid,TYPE,id,comp,q)
     
!     PUTS ID RECORD ON RANDOM OUTPUT FILES
 
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 INTEGER, INTENT(IN)                      :: mid
 INTEGER, INTENT(IN OUT)                  :: TYPE
 INTEGER, INTENT(IN)                      :: id
 INTEGER, INTENT(IN)                      :: comp
 INTEGER, INTENT(IN)                      :: q(2)
 INTEGER :: idr(50)
 
 INTEGER :: mid1(2,7)
 COMMON /output/ head(1)
 DATA mid1/2001,4HDISP, 2010,4HVELO,  &
     2011,4HACCE, 2002,4HLOAD,  &
     2003,4HSPCF, 2004,4HELFO,  &
     2005,4HSTRE/
 DATA idr /50*0/
 
 idr(1)= 50
 idr(3) = mid
 DO  i = 1,7
   IF(TYPE == mid1(2,i)) EXIT
 END DO
 20 itype = mid1(1,i)
 idr(2) = itype
 idr(5) = id*10
 idr(6) = comp
 idr(8) = q(1)
 idr(9) = q(2)
 idr(10) = 2
 CALL WRITE(FILE,idr(1),50,0)
 CALL WRITE(FILE,head(1),96,1)
 RETURN
END SUBROUTINE rand1
