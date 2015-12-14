SUBROUTINE outpt3
     
!     PUNCH UP TO 5 MATRIX DATA BLOCK ONTO DMI CARDS
 
!     CALL TO THIS MODULE IS
 
!     OUTPUT3   M1,M2,M3,M4,M5//C,N,PO/C,Y,N1=AB/C,Y,N2=CD/C,Y,N3=EF/
!                                      C,Y,N4=GH/C,Y,N5=IJ   $
 
!               PO = FORTRAN OUTPUT FILE UNIT NO. (DEFAULT = 0)
!                    .GE.0 MEANS NO LISTING OF  CARD IMAGES WILL BE MADE
!                    .LT.0 MEANS LISTING OF DMI CARD IMAGES WILL BE MADE
!                          ON FORTRAN UNIT = IABS(PO).
 
 
 
 LOGICAL :: first
 INTEGER :: in(5),subnam(2),NAME(2),trl(7),erno,param,  &
     trl1,trl2,trl3,trl4,trl5,trl6,trl7,eol,eor
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /BLANK / jo,param(2,5) /system/ nb,no,junk(6),nlpp  &
     /zzzzzz/ x(1) /zntpkx/ z(4),iz,eol,eor  &
     /phdmix/ namex(2),nam,ifo,itin,itout,ir,ic,noutpt,kpp,nlp,  &
     erno,icol,iro,xx,icard1
 EQUIVALENCE     (trl(1),trl1), (trl(2),trl2), (trl(3),trl3),  &
     (trl(4),trl4), (trl(5),trl5), (trl(6),trl6), (trl(7),trl7)
 DATA    subnam/ 4HOUTP,4HUT3 /,  in/ 101,102,103,104,105 /
 DATA    ityp  / 1 /
 
 
 lcor = korsz(x) - nb
 IF (lcor <= 0) CALL mesage (-8,lcor,subnam)
 ibuf = lcor+1
 jono = 0
 IF (jo < 0) jono = IABS(jo)
 noutpt = jono
 itin = 1
 kpp  = 2
 nlp  = nlpp
 
 DO  ii = 1,5
   trl1 = in(ii)
   CALL rdtrl (trl)
   IF (trl1 <= 0) CYCLE
   CALL fname (in(ii),NAME)
   CALL gopen (in(ii),x(ibuf),0)
   namex(1) = NAME(1)
   namex(2) = NAME(2)
   nam  = param(1,ii)
   ifo  = trl4
   itout= 0
   ir   = trl3
   ic   = trl2
   CALL phdmia
   IF (erno /= 0) GO TO 9400
   
   DO  j = 1,trl2
     CALL intpk (*900,in(ii),0,ityp,0)
     first = .false.
     icol  = j
     
     DO  i = 1,trl3
       IF (eol /= 0) GO TO 850
       CALL zntpki
       iro = iz
       xx  = z(1)
       
!     VAX MAY HAVE A FEW IMBEDED ZEROS
       
       IF (xx == 0.0) CYCLE
       IF (first) GO TO 100
       first = .true.
       CALL phdmib
       IF (erno == 0.0) THEN
         GO TO   200
       ELSE
         GO TO  9400
       END IF
       100 CALL phdmic
       IF (erno == 0.0) THEN
         GO TO   200
       ELSE
         GO TO  9400
       END IF
       200 CONTINUE
     END DO
     
     850 CALL phdmid
     IF (erno /= 0) GO TO 9400
     900 CONTINUE
   END DO
   
   ncards = icard1 + 1
   CALL page2 (-2)
   WRITE  (no,1) uim,NAME,ncards
   1 FORMAT (a29,' 4103, OUTPUT3 HAS PUNCHED MATRIX DATA BLOCK ',2A4,  &
       ' ONTO ',i5,' DMI CARDS.')
   CALL CLOSE (in(ii),1)
 END DO
 RETURN
 
!     ERROR MESSAGE
 
 9400 CALL page2 (-2)
 WRITE  (no,9450) ufm
 9450 FORMAT (a23,' 4104, ATTEMPT TO PUNCH MORE THAN 99999 DMI CARDS ',  &
     'FOR A SINGLE MATRIX.')
 CALL mesage (-61,0,0)
 RETURN
 
END SUBROUTINE outpt3
