SUBROUTINE flbprt (iuset,ieqex,ibuf)
     
!     HYDROELEASTIC USET OUTPUT
 
!     PRINTS DOF VS. DISP SETS IF DIAG 32 IS ON.
!     PRINTS DISP SETS VS. DOF IF DIAG 33 IS ON.
 
 
 INTEGER, INTENT(IN)                      :: iuset
 INTEGER, INTENT(IN)                      :: ieqex
 INTEGER, INTENT(IN)                      :: ibuf
 EXTERNAL        andf
 INTEGER :: z        ,sysbuf   ,eqexin   ,d32      ,d33      ,  &
     FILE     ,NAME(2)  ,msk(17)  ,zgrd(10) ,title(3,9)  &
     ,               zdof(10) ,two      ,um       ,uo       ,ur       ,  &
     usg      ,usb      ,ul       ,ua       ,uf       ,  &
     us       ,un       ,ug       ,ux       ,uy       ,  &
     ufr      ,uz       ,uab      ,ui       ,dash     ,  &
     astric   ,BLANK    ,sbit(17) ,expnt    ,upbit(17), andf
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm
 COMMON /flbfil/ dum1(8)  ,eqexin
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf   ,nout     ,dum2(6)  ,nlpp     ,  &
     dum3     ,npage    ,line
 COMMON /two   / two(32)
 COMMON /bitpos/ um       ,uo       ,ur       ,usg      ,  &
     usb      ,ul       ,ua       ,uf       ,  &
     us       ,un       ,ug       ,ue       ,  &
     up       ,une      ,ufe      ,ud       ,  &
     ups      ,usa      ,uk       ,upa      ,  &
     u21      ,u22      ,u23      ,ux       ,  &
     uy       ,ufr      ,uz       ,uab      , ui
 DATA    NAME  / 4HFLBP   , 4HRT    /
 DATA    title / 4H       , 4H      , 4H mpc  ,  &
     4H       , 4H      , 4H spc  , 4H       , 4H      , 4HOMIT  ,  &
     4H       , 4HANAL  , 4HYSIS  , 4H       , 4HPERM  , 4H spc  ,  &
     4H       , 4HBDRY  , 4H spc  , 4H   s   , 4HTRUC  , 4HTURE  ,  &
     4H       , 4H   f  , 4HLUID  , 4HFREE   , 4H sur  , 4HFACE  /
 DATA    BLANK / 1H  /    , dash    / 1H- /   , astric / 1H* /
 
 
!     DETERMINE IF ANY OUTPUT IS REQUESTED
 
 CALL sswtch (32,d32)
 CALL sswtch (33,d33)
 IF (d32 == 0 .AND. d33 == 0) RETURN
 
!     READ EQEXIN INTO CORE
 
 FILE = eqexin
 CALL OPEN (*1001,eqexin,z(ibuf),0)
 CALL fwdrec (*1002,eqexin)
 CALL fwdrec (*1002,eqexin)
 nz = ibuf - ieqex
 CALL READ (*1002,*10,eqexin,z(ieqex),nz,1,neqex)
 GO TO 1008
 10 CALL CLOSE (eqexin,1)
 
!     SORT ON INTERNAL ID
 
 CALL sort (0,0,2,2,z(ieqex),neqex)
 
!     SET UP USET MASKS FOR DOF VS. DISP SET PRINTOUT
 
 IF (d32 == 0) GO TO 100
 msk( 1) = two(usb)
 msk( 2) = two(usg)
 msk( 3) = two( ul)
 msk( 4) = two( ua)
 msk( 5) = two( uf)
 msk( 6) = two( un)
 msk( 7) = two( ug)
 msk( 8) = two( ur)
 msk( 9) = two( uo)
 msk(10) = two( us)
 msk(11) = two( um)
 msk(12) = two( ux)
 msk(13) = two( uy)
 msk(14) = two(ufr)
 msk(15) = two( uz)
 msk(16) = two(uab)
 msk(17) = two( ui)
 
 DO  i = 1,17
   sbit(i) = 0
 END DO
 
!     PASS THROUGH EQEXIN TABLE AND DETERMINE NUMBER OF DOF FOR EACH
!     POINT
 
 juset = iuset - 1
 line  = nlpp
 inpnt = 0
 DO  k = 1,neqex,2
   itype = MOD(z(ieqex+k),10)
   ndof  = 6
   IF (itype == 2) ndof = 1
   
!     FOR EACH DOF - GET USET ENTRY AND TEST VARIOUS MACK BITS
   
   DO  kk = 1,ndof
     juset = juset + 1
     iu    = z(juset)
     inpnt = inpnt + 1
     expnt = z(ieqex+k-1)
     idof  = kk
     IF (ndof == 1) idof = 0
     DO  ibit = 1,17
       IF (andf(msk(ibit),iu) /= 0) GO TO 25
       upbit(ibit) = BLANK
       CYCLE
       25 upbit(ibit) = astric
       sbit (ibit) = sbit(ibit) + 1
     END DO
     
!     PRINT LINE OF OUTPUT
     
     line = line + 1
     IF (line <= nlpp) GO TO 40
     CALL page1
     WRITE (nout,2000)
     line = 1
     40 WRITE (nout,2010) inpnt,expnt,dash,idof,upbit
   END DO
 END DO
 
!     PRINT COLUMN TOTALS
 
 WRITE (nout,2020) sbit
 
!     SET UP MASKS FOR DISP SET VS. DOF PRINTOUT
 
 100 IF (d33 == 0) RETURN
 msk( 1) = two( um)
 msk( 2) = two( us)
 msk( 3) = two( uo)
 msk( 4) = two( ua)
 msk( 5) = two(usg)
 msk( 6) = two(usb)
 msk( 7) = two( ux)
 msk( 8) = two( uy)
 msk( 9) = two(ufr)
 
!     PASS THROUGH EQEXIN TABLE ONCE FOR EACH DISP SET TO BE PRINTED
 
 DO  imk = 1,9
   inum  = -9
   icol  =  0
   line  = nlpp
   juset = iuset - 1
   DO  k = 1,neqex,2
     itype = MOD(z(ieqex+k),10)
     ndof  = 6
     IF (itype == 2) ndof = 1
     
!     FOR EACH DOF - TEST IF IT IS IN DESIRED SET FOR THIS PASS
     
     expnt = z(ieqex+k-1)
     DO  kk = 1,ndof
       juset = juset + 1
       IF (andf(z(juset),msk(imk)) == 0) CYCLE
       idof = kk
       IF (ndof == 1) idof = 0
       icol = icol + 1
       zgrd(icol) = expnt
       zdof(icol) = idof
       IF (icol < 10) CYCLE
       
!     WE HAVE ACUMULATED 10 POINTS - PRINT THEM
       
       icol = 0
       line = line + 1
       IF (line <= nlpp) GO TO 110
       CALL page1
       WRITE (nout,2030) (title(i,imk),i=1,3)
       line = 1
       110 inum = inum + 10
       WRITE (nout,2040) inum,(zgrd(i),zdof(i),i=1,10)
       
     END DO
   END DO
   
!     PRINT ANY REMAINING ENTRIES
   
   IF (icol == 0) CYCLE
   line = line + 1
   IF (line <= nlpp) GO TO 140
   CALL page1
   WRITE (nout,2030) (title(i,imk),i=1,3)
   line = 1
   140 inum = inum + 10
   WRITE (nout,2040) inum,(zgrd(i),zdof(i),i=1,icol)
   
 END DO
 
!     PRINT OUT COMPLETE
 
 RETURN
 
!     ERROR CONDITIONS - PRINT NON-FATAL MESSAGE
 
 1001 n = 1
 GO TO 1100
 1002 n = 2
 GO TO 1100
 1008 WRITE (nout,2050) uwm
 RETURN
 1100 CALL mesage (n,FILE,NAME)
 RETURN
 
!     FORMAT STATEMENTS
 
 2000 FORMAT (//12X,'INT DOF  EXT GP. DOF   SB   SG    L    A    F   ',  &
     'N    G    R    O    S    M    X    Y   FR    Z   AB    I', /1X,131(1H-))
 2010 FORMAT (10X,i8,1X,i8,1X,a1,i2,1X,17(4X,a1))
 2020 FORMAT (1H0,31H-- c o l u m n   t o t a l s -- ,17I5)
 2030 FORMAT (45X,3A4,17H displacement set, //16X,3H-1-,8X,3H-2-,8X,  &
     3H-3-,8X,3H-4-,8X,3H-5-,8X,3H-6-,8X,3H-7-,8X,3H-8-,8X, 3H-9-,7X,4H-10- ,/1H )
 2040 FORMAT (1H ,i6,1H=,10(1X,i8,1H-,i1))
 2050 FORMAT (a25,' 8011, INSUFFICIENT CORE TO HOLD CONTENTS OF EQEXIN',  &
     ' DATA BLOCK', /31X,'HYDROELASTIC USET PRINTOUT TERMINATED.')
END SUBROUTINE flbprt
