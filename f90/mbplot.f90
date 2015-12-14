SUBROUTINE mbplot (nw1,nd1,nwn,nc21,nc2n,nc1,ncn,ndn)
     
!     SUBROUTINE TO PRINT A REPRESENTATION OF PLANFORM BEING CONSIDERED
 
 
 INTEGER, INTENT(IN)                      :: nw1(1)
 INTEGER, INTENT(IN)                      :: nd1(1)
 INTEGER, INTENT(IN OUT)                  :: nwn(1)
 INTEGER, INTENT(IN)                      :: nc21(1)
 INTEGER, INTENT(IN)                      :: nc2n(1)
 INTEGER, INTENT(IN)                      :: nc1(1)
 INTEGER, INTENT(IN OUT)                  :: ncn(1)
 INTEGER, INTENT(IN OUT)                  :: ndn(1)
 REAL :: mach
 DIMENSION  pl(50)
 COMMON /system/ sys,n6
 COMMON /mboxc / njj,crank1,crank2,cntrl1,cntrl2,nbox,npts0,npts1,  &
     npts2,asym,gc,cr,mach,beta,ek,ekbar,ekm,boxl,boxw,  &
     boxa ,ncb,nsb,nsbd,ntote,kc,kc1,kc2,kct,kc1t,kc2t
 DATA    BLANK , dia, wg , fp , tp , wk  / 1H    , 1H., 1HS, 1H1,1H2 , 1H, /
 
 nsbm  = MAX0(nsb,nsbd)
 ncbmx = MAX0(ncb,5   )
 WRITE  (n6,200) mach,boxw,boxl
 200  FORMAT (1H1,29X,'GRAPHIC DISPLAY OF REGIONS ON MAIN SEMISPAN',  &
     /10X,11HMACH NUMBER ,f8.3,11X,9HBOX width ,f11.6 ,10X,  &
     10HBOX length ,f11.6, //)
 DO  i = 1,ncbmx
   DO  j = 1,nsbm
     pl(j) = BLANK
     IF (j > nsb   ) GO TO 1500
     IF (i >= nw1(j)) GO TO 1100
     IF (i < nd1(j)) CYCLE
     pl(j) = dia
     CYCLE
     1100 IF (i > nwn(j)) GO TO 1300
     IF (i >= nc21(j) .AND. i <= nc2n(j)) GO TO 1150
     IF (i >= nc1(j)  .AND. i <= ncn(j) ) GO TO 1200
     pl(j) = wg
     CYCLE
     1150 pl(j) = tp
     CYCLE
     1200 pl(j) = fp
     CYCLE
     1300 IF (i > ndn(j)) CYCLE
     pl(j) = wk
     CYCLE
     1500 IF ((i >= nd1(j) .AND. i <= ndn(j)) .OR. (i >= nc1(j) .AND.  &
         i <= ncn(j))) pl(j) = dia
     CYCLE
   END DO
   
   WRITE  (n6,2000) (pl(j),j=1,nsbm)
   2000 FORMAT (30X,50A1)
   
   IF (i > 5) CYCLE
   SELECT CASE ( i )
     CASE (    1)
       GO TO 2100
     CASE (    2)
       GO TO 2300
     CASE (    3)
       GO TO 2500
     CASE (    4)
       GO TO 2700
     CASE (    5)
       GO TO 2900
   END SELECT
   2100 WRITE  (n6,2200)
   2200 FORMAT (1H+,84X,9HS    main )
   CYCLE
   2300 WRITE  (n6,2400)
   2400 FORMAT (1H+,84X,11H1    cntrl1 )
   CYCLE
   2500 WRITE  (n6,2600)
   2600 FORMAT (1H+,84X,11H2    cntrl2 )
   CYCLE
   2700 WRITE  (n6,2800)
   2800 FORMAT (1H+,84X,14H.    diaphragm )
   CYCLE
   2900 WRITE  (n6,3000)
   3000 FORMAT (1H+,84X,9H,    wake )
 END DO
 RETURN
END SUBROUTINE mbplot
