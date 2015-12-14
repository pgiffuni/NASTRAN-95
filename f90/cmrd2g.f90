SUBROUTINE cmrd2g
     
!     THIS SUBROUTINE CREATES THE REDUCED SUBSTRUCTURE NEW TABLE ITEMS
!     FOR THE CMRD2 MODULE.
 
!     INPUT DATA
!     GINO  - EQST  - TEMPORARY SUBSTRUCTURE EQUIVALENCE TABLE FOR
!                     SUBSTRUCTURE BEING REDUCED
 
!     OUTPUT DATA
!     SOF  - EQSS   - SUBSTRUCTURE EQUIVALENCE TABLE FOR REDUCED
!                     SUBSTRUCTURE
!            BGSS   - BASIC GRID POINT DEFINITION TABLE FOR REDUCED
!                     SUBSTRUCTURE
!            LODS   - LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!            LOAP   - APPENDED LOAD SET DATA FOR REDUCED SUBSTRUCTURE
!            PLTS   - PLOT SET DATA FOR REDUCED SUBSTRUCTURE
!            CSTM   - COORDINATE SYSTEM TRANSFORMATION DATA FOR REDUCED
!                     SUBSTRUCTURE
 
!     PARAMETERS
!     INPUT- DRY    - MODULE OPERATION FLAG
!            POPT   - LOADS OPERATION FLAG
!            GBUF1  - GINO BUFFER
!            INFILE - INPUT FILE NUMBERS
!            KORLEN - LENGTH OF OPEN CORE
!            KORBGN - BEGINNING ADDRESS OF OPEN CORE
!            OLDNAM - NAME OF SUBSTRUCTURE BEING REDUCED
!            NEWNAM - NAME OF REDUCED SUBSTRUCTURE
!            FREBDY - FREEBODY OPR
!            FREBDY - FREEBODY OPTIONS FLAG
!            IO     - OUTPUT OPTIONS FLAG
!            MODPTS - NUMBER OF MODAL POINTS
 
 EXTERNAL        rshift,andf
 LOGICAL :: ponly
 INTEGER :: dry,popt,gbuf1,oldnam,z,andf,rshift,eqst,sofeog
 DIMENSION       modnam(2),lstbit(32),itrlr(7),itmlst(3),itmnam(2), rz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / idum1,dry,popt,gbuf1,idum2(5),infile(11),  &
     idum3(17),korlen,korbgn,oldnam(2),newnam(2),  &
     idum4(4),io,idum5(3),modpts,idum9,ponly
 COMMON /zzzzzz/ z(1)
 COMMON /system/ idum6,iprntr,idum7(6),nlpp,idum8(2),line
 EQUIVALENCE     (eqst,infile(5)),(rz(1),z(1))
 DATA    modnam/ 4HCMRD,4H2G  /
 DATA    papp  , lods,loap    /4HPAPP,4HLODS,4HLOAP   /
 DATA    itmlst/ 4HEQSS,4HBGSS,4HLAMS   /
 DATA    sofeog/ 4H$eog/, nhplts,nhcstm /4HPLTS,4HCSTM/
 
!     CHECK FOR LOADS ONLY
 
 IF (ponly) GO TO 55
 
!     PROCESS EQSS, BGSS DATA
 
 IF (dry == -2) RETURN
 itrlr(1) = eqst
 CALL rdtrl (itrlr)
 ifile = eqst
 IF (itrlr(1) < 0) GO TO 210
 CALL gopen (eqst,z(gbuf1),0)
 itest = 3
 item  = itmlst(1)
 itmnam(1) = newnam(1)
 itmnam(2) = newnam(2)
 CALL sfetch (newnam,item,2,itest)
 IF (itest /= 3) GO TO 250
 newpts = modpts
 
!     PROCESS EQSS GROUP 0 DATA
 
 IF (korbgn+itrlr(2)+2 >= korlen) GO TO 230
 CALL READ (*215,*220,eqst,z(korbgn),itrlr(2),1,nwdsrd)
 ncsubs = z(korbgn+2)
 z(korbgn+2) = z(korbgn+2) + 1
 z(korbgn+3) = z(korbgn+3) + newpts
 newcs  = itrlr(2)
 z(korbgn+newcs  ) = newnam(1)
 z(korbgn+newcs+1) = newnam(2)
 newcs  = itrlr(2) + 2
 CALL suwrt (z(korbgn),newcs,2)
 
!     PROCESS REMAINING EQSS GROUPS
 
 nwds = korlen - korbgn
 DO  i = 1,ncsubs
   CALL READ (*215,*10,eqst,z(korbgn),nwds,1,nwdsrd)
   GO TO 230
   10 IF (korbgn+1+nwdsrd >= korlen) GO TO 230
   CALL suwrt (z(korbgn),nwdsrd,2)
 END DO
 
!     PROCESS MODAL POINTS
 
 IF (korbgn+3*newpts >= korlen) GO TO 230
 DO  i = 1,newpts
   kore = 3*(i-1)
   z(korbgn+kore  ) = 100 + i
   z(korbgn+kore+1) = itrlr(4)/2 + i
   z(korbgn+kore+2) = 1
 END DO
 nwdsrd = 3*newpts
 CALL suwrt (z(korbgn),nwdsrd,2)
 
!     PROCESS EQSS SIL DATA
 
 IF (korbgn+itrlr(4)+2*newpts >= korlen) GO TO 230
 CALL READ (*215,*220,eqst,z(korbgn),itrlr(4),1,nwdsrd)
 nwdsrd = itrlr(4) - 1
 icode  = z(korbgn+nwdsrd)
 CALL decode (icode,lstbit,nwdsd)
 lstsil = z(korbgn+nwdsrd-1) + nwdsd - 1
 DO  i = 1,newpts
   kore = itrlr(4) + 2*(i-1)
   z(korbgn+kore  ) = lstsil + i
   z(korbgn+kore+1) = 1
 END DO
 nwdsrd = itrlr(4) + 2*newpts
 CALL suwrt (z(korbgn),nwdsrd,2)
 CALL suwrt (z(korbgn),0,3)
 
!     PROCESS BGSS DATA
 
 IF (korbgn+itrlr(5)+4*newpts >= korlen) GO TO 230
 item  = itmlst(2)
 itest = 3
 CALL sfetch (newnam,item,2,itest)
 IF (itest /= 3) GO TO 250
 CALL READ (*215,*220,eqst,z(korbgn),3,1,nwdsrd)
 z(korbgn  ) = newnam(1)
 z(korbgn+1) = newnam(2)
 z(korbgn+2) = z(korbgn+2) + newpts
 locbgs = korbgn
 CALL suwrt (z(korbgn),3,2)
 CALL READ (*215,*220,eqst,z(korbgn),itrlr(5),1,nwdsrd)
 DO  i = 1,newpts
   kore = itrlr(5) + 4*(i-1)
   z(korbgn+kore   ) = -1
   rz(korbgn+kore+1) = 0.0
   rz(korbgn+kore+2) = 0.0
   rz(korbgn+kore+3) = 0.0
 END DO
 nwdsrd = itrlr(5) + 4*newpts
 CALL suwrt (z(korbgn),nwdsrd,2)
 CALL suwrt (z(korbgn),0,3)
 korbgn = korbgn + itrlr(5)
 
!     PROCESS LODS, LOAP ITEM
 
 55 item = lods
 IF (popt == papp) item = loap
 itest = 3
 CALL sfetch (oldnam,item,1,itest)
 IF (itest == 3) GO TO 60
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 IF (korbgn+nwdsrd >= korlen) GO TO 230
 z(korbgn  ) = newnam(1)
 z(korbgn+1) = newnam(2)
 z(korbgn+3) = z(korbgn+3) + 1
 z(korbgn+nwdsrd  ) = newnam(1)
 z(korbgn+nwdsrd+1) = newnam(2)
 z(korbgn+nwdsrd+2) = sofeog
 iwds = nwdsrd + 3
 CALL suread (z(korbgn+iwds),-2,nwdsrd,itest)
 IF (korbgn+iwds+nwdsrd+2 >= korlen) GO TO 230
 z(korbgn+iwds+nwdsrd  ) = 0
 z(korbgn+iwds+nwdsrd+1) = sofeog
 iwds  = iwds + nwdsrd + 2
 itest = 3
 CALL sfetch (newnam,item,2,itest)
 IF (itest /= 3) GO TO 250
 CALL suwrt (z(korbgn),iwds,3)
 IF (ponly) GO TO 130
 
!     PROCESS PLTS ITEM
 
 60 CALL sfetch (oldnam,nhplts,1,itest)
 IF (itest == 3) GO TO 70
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 z(korbgn  ) = newnam(1)
 z(korbgn+1) = newnam(2)
 itest = 3
 CALL sfetch (newnam,nhplts,2,itest)
 IF (itest /= 3) GO TO 250
 CALL suwrt (z(korbgn),nwdsrd,itest)
 
!     PROCESS CSTM ITEM
 
 70 CALL sfetch (oldnam,nhcstm,1,itest)
 IF (itest == 3) GO TO 130
 CALL suread (z(korbgn),-2,nwdsrd,itest)
 IF (korbgn+2*nwdsrd >= korlen) GO TO 230
 z(korbgn  ) = newnam(1)
 z(korbgn+1) = newnam(2)
 kore = nwdsrd - 4
 CALL sort (0,0,14,1,z(korbgn+3),kore)
 kore = kore/14
 IF (korbgn+2*nwdsrd+kore >= korlen) GO TO 230
 DO  i = 1,kore
   z(korbgn+nwdsrd+i-1) = 0
 END DO
 nbgss = itrlr(5)/4
 loop100:  DO  i = 1,nbgss
   k  = 4*(i-1)
   IF (z(locbgs+k) <= 0) CYCLE loop100
   DO  j = 1,kore
     loc = 14*(j-1)
     IF (z(korbgn+3+loc) /= z(locbgs+k)) CYCLE
     z(korbgn+nwdsrd+j-1) = 1
     CYCLE loop100
   END DO
 END DO loop100
 locnew = 0
 DO  i = 1,kore
   IF (z(korbgn+nwdsrd+i-1) == 0) CYCLE
   locold = 14*(i-1)
   DO  j = 1,14
     z(korbgn+nwdsrd+kore+locnew+j-1) = z(korbgn+3+locold+j-1)
   END DO
   locnew = locnew + 14
 END DO
 IF(locnew == 0) GO TO 130
 itest = 3
 CALL sfetch (newnam,nhcstm,2,itest)
 CALL suwrt (newnam,2,2)
 CALL suwrt (z(korbgn+nwdsrd+kore),locnew,2)
 CALL suwrt (z(korbgn),0,3)
 
!     OUTPUT EQSS ITEM
 
 130 CALL CLOSE (eqst,1)
 IF (andf(rshift(io,4),1) /= 1) GO TO 150
 CALL sfetch (newnam,itmlst(1),1,itest)
 IF (itest /= 1) GO TO 250
 CALL suread (z(korbgn),4,nwdsrd,itest)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 loc    = korbgn + nwdsrd
 ncsubs = ncsubs + 1
 DO  i = 1,ncsubs
   CALL suread (z(loc),-1,nwdsrd,itest)
   namloc = korbgn + 2*(i-1)
   CALL cmiwrt (1,newnam,z(namloc),loc,nwdsrd,z,z)
 END DO
 CALL suread (z(loc),-1,nwdsrd,itest)
 IF (loc+nwdsrd >= korlen) GO TO 230
 CALL cmiwrt (8,newnam,0,loc,nwdsrd,z,z)
 
!     OUTPUT BGSS ITEM
 
 150 IF (andf(rshift(io,5),1) /= 1) GO TO 160
 CALL sfetch (newnam,itmlst(2),1,itest)
 IF (itest /= 1) GO TO 250
 ngrp = 1
 CALL sjump (ngrp)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL cmiwrt (2,newnam,newnam,korbgn,nwdsrd,z,z)
 
!     OUTPUT CSTM ITEM
 
 160 IF (andf(rshift(io,6),1) /= 1) GO TO 170
 CALL sfetch (newnam,nhcstm,1,itest)
 IF (itest == 3) GO TO 170
 ngrp = 1
 CALL sjump (ngrp)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL cmiwrt (3,newnam,newnam,korbgn,nwdsrd,z,z)
 
!     OUTPUT PLTS ITEM
 
 170 IF (andf(rshift(io,7),1) /= 1) GO TO 180
 CALL sfetch (newnam,nhplts,1,itest)
 IF (itest == 3) GO TO 180
 CALL suread (z(korbgn),3,nwdsrd,itest)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL cmiwrt (4,newnam,newnam,korbgn,nwdsrd,z,z)
 
!     OUTPUT LODS ITEM
 
 180 IF (andf(rshift(io,8),1) /= 1) GO TO 200
 CALL sfetch (newnam,  lods,1,itest)
 IF (itest == 3) GO TO 200
 CALL suread (z(korbgn),4,nwdsrd,itest)
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 loc   = korbgn + nwdsrd
 itype = 5
 IF (item == loap) itype = 7
 DO  i = 1,ncsubs
   namloc = korbgn + 2*(i-1)
   CALL suread (z(loc),-1,nwdsrd,itest)
   CALL cmiwrt (itype,newnam,z(namloc),loc,nwdsrd,z,z)
   itype = 6
 END DO
 
!     OUTPUT MODAL DOF SUMMARY
 
 200 IF (andf(rshift(io,9),1) /= 1) GO TO 209
 item = itmlst(3)
 itmnam(1) = oldnam(1)
 itmnam(2) = oldnam(2)
 CALL sfetch (oldnam,item,1,itest)
 IF (itest /= 1) GO TO 250
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 CALL page1
 WRITE (iprntr,901) newnam
 line   = line + 10
 nofreq = z(korbgn+3)
 lamloc = korbgn
 moduse = lamloc + 7*nofreq + 1
 CALL suread (z(korbgn),-2,nwdsrd,itest)
 IF (korbgn+nwdsrd >= korlen) GO TO 230
 item   = itmlst(1)
 itmnam(1) = newnam(1)
 itmnam(2) = newnam(2)
 CALL sfetch (newnam,item,1,itest)
 IF (itest /= 1) GO TO 250
 korbgn = korbgn + moduse + nofreq
 IF (korbgn >= korlen) GO TO 230
 CALL suread (z(korbgn),-1,nwdsrd,itest)
 DO  i = 1,ncsubs
   CALL suread (z(korbgn),-1,nwdsrd,itest)
   IF (korbgn+nwdsrd >= korlen) GO TO 230
 END DO
 loceqs = korbgn
 ipid   = 2*z(korbgn+1)
 korbgn = korbgn + nwdsrd
 IF (korbgn+ipid >= korlen) GO TO 230
 CALL suread (z(korbgn),ipid,nwdsrd,itest)
 ips    = z(korbgn+ipid-2)
 index1 = -3
 DO  i = 1,nofreq
   IF (line <= nlpp) GO TO 204
   CALL page1
   WRITE (iprntr,901) newnam
   line = line + 10
   204 CONTINUE
   IF (z(moduse+i-1) > 1) GO TO 206
   index1 = index1 + 3
   mode   = 7*(i-1)
   WRITE (iprntr,902) z(lamloc+mode),rz(lamloc+mode+4),z(moduse+i-1),  &
       z(loceqs+index1),ips
   ips = ips + 1
   GO TO 208
   206 mode = 7*(i-1)
   WRITE (iprntr,902) z(lamloc+mode),rz(lamloc+mode+4),z(moduse+i-1)
   208 line = line + 1
 END DO
 209 CONTINUE
 RETURN
 
!     PROCESS SYSTEM FATAL ERRORS
 
 210 imsg = -1
 GO TO 240
 215 imsg = -2
 GO TO 240
 220 imsg = -3
 GO TO 240
 230 imsg = -8
 ifile = 0
 240 CALL sofcls
 CALL mesage (imsg,ifile,modnam)
 RETURN
 
!     PROCESS MODULE FATAL ERRORS
 
 250 SELECT CASE ( itest )
   CASE (    1)
     GO TO 260
   CASE (    2)
     GO TO 260
   CASE (    3)
     GO TO 260
   CASE (    4)
     GO TO 270
   CASE (    5)
     GO TO 280
   CASE (    6)
     GO TO 280
 END SELECT
 260 WRITE (iprntr,900) ufm,modnam,item,itmnam
 dry = -2
 RETURN
 
 270 imsg = -2
 GO TO 290
 280 imsg = -3
 290 CALL smsg (imsg,item,itmnam)
 RETURN
 
 900 FORMAT (a23,' 6211, MODULE ',2A4,' - ITEM ',a4,  &
     ' OF SUBSTRUCTURE ',2A4,' HAS ALREADY BEEN WRITTEN.')
 901 FORMAT (//36X,43HMODAL dof summary for reduces substructure ,2A4,  &
     //30X,41HUSAGE codes are 1 - included in modal set,  &
     /46X,50H2 - excluded from modal set because of non-partici,  &
     6HPATION,/46X,41H3 - excluded from modal set because of ra,  &
     11HNGE OR nmax, //40X,4HMODE,22X,15HUSAGE      grid, /39X,  &
     6HNUMBER,8X,6HCYCLES,8X,26HCODE    point id       sil,/)
 902 FORMAT (39X,i5,5X,1P,e13.6,6X,i1,6X,i8,4X,i6)
 
END SUBROUTINE cmrd2g
