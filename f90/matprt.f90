SUBROUTINE matprt (*,*,a,option,column)
     
!     MATPRT AND PRTMAT ARE CALLED ONLY BY INTPRT
 
 REAL, INTENT(IN OUT)                     :: a(1)
 INTEGER, INTENT(IN OUT)                  :: option
 INTEGER, INTENT(OUT)                     :: column(1)
 
 INTEGER :: FILE,TYPE,bufsiz,count, utype,ui,uj,uinc,rsp,rdp,csp,cdp,rew
 COMMON /system/ bufsiz,mo,skp1(6),maxlin,skp2(2),count
 COMMON /unpakx/ utype,ui,uj,uinc
 COMMON /xxmprt/ mcb(7)
 
!     MCB = MATRIX CONTROL BLOCK.
!     A   = ARRAY OF BUFSIZ + I (REAL) OR 2I (COMPLEX) LOCATIONS.
!     OPTION IS AS DESCRIBED IN -VECPRT-.
!     RETURN 1 ... PRINT MATRIX TITLE + COLUMN IDENTIFIER.
!     RETURN 2 ... PRINT COLUMN IDENTIFIER ONLY.
!                  (PRTMAT = RETURN ENTRY POINT)
!     COLUMN = CURRENT COLUMN NUMBER
 
 EQUIVALENCE (FILE,mcb(1)), (j,mcb(2)), (i,mcb(3)), (TYPE,mcb(5))
 DATA        rsp,rdp,csp,cdp,rew,inprew / 1,2,3,4,1,0 /
 
 IF (i <= 0 .OR. j <= 0)  GO TO 150
 utype = TYPE
 IF (TYPE == rdp) utype = rsp
 IF (TYPE == cdp) utype = csp
 ui   = 1
 uj   = i
 uinc = 1
 CALL gopen (FILE,a,inprew)
 count = maxlin
 
 column(1) = 0
 110 column(1) = column(1) + 1
 CALL unpack (*140,FILE,a(bufsiz+1))
 CALL vecprt (*120,*130,utype,i,a(bufsiz+1),option)
 GO TO 140
 120 RETURN 1
 130 RETURN 2
 
 
 ENTRY prtmat (*,*,column)
!     =========================
 
 CALL prtvec (*120,*130)
 140 IF (column(1) /= j) GO TO 110
 
 CALL CLOSE (FILE,rew)
 
 150 RETURN
END SUBROUTINE matprt
