SUBROUTINE ddr1
     
!     DYNAMIC DATA RECOVERY PART1
 
!     INPUTS  2 UHV,PHIDH
 
!     OUTPUTS 1 UDV
 
!     SCRATCHES  1
 
 INTEGER :: uhv,phidh,udv,scr1
 DATA uhv,phidh,udv,scr1/101,102,201,301/
 
!     TRANSFPRM TO MODAL DISPLACEMENTS
 
 CALL ssg2b(phidh,uhv,0,udv,0,1,1,scr1)
 RETURN
END SUBROUTINE ddr1
