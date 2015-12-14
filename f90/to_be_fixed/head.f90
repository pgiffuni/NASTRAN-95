SUBROUTINE head (dtyp,pltp,mtyp,idat)
     
 
 INTEGER, INTENT(IN OUT)                  :: dtyp
 INTEGER, INTENT(IN OUT)                  :: pltp
 INTEGER, INTENT(IN OUT)                  :: mtyp
 INTEGER, INTENT(OUT)                     :: idat(17)
 INTEGER :: maxdef(3), undef(4),ptyp(2,5),  &
     subc(2),mtypf(2,3),phas(3),fpltit,pltitl
 REAL :: nt1(5),nt2(4),nt3(3),cscale,x,x0
 COMMON /output/ title(32,3)
 COMMON /pltdat/ skpplt(2),xymin(2),xymax(2),axymax(13),cscale,  &
     skpa(3),cntx,cnty
 COMMON /xxparm/ iskp(215),fpltit,pltitl(17)
 
 DATA    undef / 4HUNDE, 4HFORM, 4HED s, 4HHAPE /
! ... NUMBER CHAR+2 FOR STATIC - CMODAL ... NOTE, 1 BLANK AT START...  &
 ,       nt1   / 8., 7., 8., 7., 8.  /
! ... NUMBER CHAR+1 FOR DEFO - ACCEL ...  &
 ,       nt2   , ptyp / 7., 9., 7., 7.  &
     ,               4HDEFO,2HR.  , 4HVELO,4HCITY, 4HACCE,2HL.  &
     ,               4HSTRE,2HSS  , 4HSTRA,2HIN  / ,       subc  / 4HSUBC,4HASE  /
! ... NUMBER CHAR+1 FOR FREQ, EIGENV., TIME  ... IDENTIFY BY MTYP ...  &
 ,       nt3   / 6., 8., 5. /  &
     ,       mtypf / 4HFREQ,4H.   , 4HEIGE,4HNV. , 4HTIME,1H    /
! ... NUMBER OF SPACES BETWEEN IDENTIFIERS ...  &
 ,       delx  / 3.0 / ,       maxdef/ 4HMAX-,4HDEF.,2H =   /  &
     ,       phas  / 4H pha,4HSE  ,1H     /
 
 xymin(1) = 0.0
 xymin(2) = 0.0
 xymax(1) = axymax(1)
 xymax(2) = axymax(2)
 CALL PRINT (0,0,0,0,0,-1)
 IF (mtyp < 0) GO TO 30
 
!     LEFT-MOST CHARACTER MAY NOT BE COMPETELY DRAWN IF FRACTION OF
!     CSCALE IS IS LESS THAN 0.5. SO MOVE OVER A SMALL SPACE OF X0
 
 j  = IFIX(cscale)
 x0 = cscale - FLOAT(j)
 IF (x0 > 0.5) x0 = 0.0
 
!     PRINT THE TITLE, SUBTITLE AND LABEL
 
 CALL PRINT (x0,3.0*cnty,1,title(1,1),17,0)
 CALL PRINT (x0,2.0*cnty,1,title(1,2),16,0)
 CALL PRINT (x0,cnty,1,title(1,3),17,0)
 
 x = 25. - 5.*(cscale-1.)
 IF (dtyp == 0) GO TO 10
 x = 40.
 IF (idat(1) <= 8) GO TO 10
 x = 45.
 IF (idat(1) >= 12) x = 52.
 IF (idat(1) >= 15) x = 59.
 10 CONTINUE
 IF (fpltit /= 0) CALL PRINT (x*cntx,0.,1,pltitl,17,0)
 
!     BOTTOM LINE IDENTIFIES PLOT
 
 IF (dtyp /= 0) GO TO 20
 
!     UNDEFORMED SHAPE
 
 CALL PRINT (cntx+x0,0.,1,undef,4,0)
 GO TO 40
 
!     DEFORMED SHAPE
 
 20 CALL PRINT (cntx+x0,0.,1,idat(3),2,0)
 x = nt1(dtyp)
 CALL PRINT (x*cntx+x0,0.,1,ptyp(1,pltp),2,0)
 x = x + nt2(pltp)
 CALL PRINT (x*cntx+x0,0.,1,subc,2,0)
 x = x + 8.
 n = -1
 CALL typint (x*cntx+x0,0.,1,idat(7),n,0)
 x = x + FLOAT(n) + delx
 
!     LOAD I  OR  MODE I
 
 CALL PRINT (x*cntx+x0,0.,1,idat(9),1,0)
 x = x + 5.
 n = -1
 CALL typint (x*cntx+x0,0.,1,idat(8),n,0)
 
!     FREQUENCY, EIGENVALUE, OR TIME
 
 IF (idat(1) <= 8) GO TO 40
 x = FLOAT(IFIX(x+delx+0.1) + n)
 CALL PRINT (x*cntx+x0,0.,1,mtypf(1,mtyp),2,0)
 x = x + nt3(mtyp)
 CALL typflt (x*cntx+x0,0.,1,idat(10),-8,0)
 
!     MAGNITUDE  OR  PHASE LAG
 
 IF (idat(1) <= 12) GO TO 40
 x = x + 7.0 + delx
 IF (idat(14) /= phas(1)) GO TO 25
 idat(15) = phas(2)
 idat(16) = phas(3)
 25 CALL PRINT (x*cntx+x0,0.,1,idat(14),3,0)
 
 IF (idat(1) <= 15) GO TO 40
 x = x + 7.0
 CALL typflt (x*cntx+x0,0.,1,idat(17),-6,0)
 GO TO 40
 
!     PRINT THE MAXIMUM DEFORMATION AT THE TOP
 
 30 CALL PRINT (20.*cntx,xymax(2),1,maxdef,3,0)
 CALL typflt (31.*cntx,xymax(2),1,idat(1),-10,0)
 
 
 40 CALL PRINT (0,0,0,0,0,1)
 RETURN
END SUBROUTINE head
