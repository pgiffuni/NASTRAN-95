SUBROUTINE setlvl (newnm,numb,oldnms,itest,ibit)
     
!     CREATES A NEW SUBSTRUCTURE NEWNM WHERE
!     - NEWNM IS AN INDEPENDENT SUBSTRUCTURE IF NUMB = 0
!     - NEWNM IS REDUCED FROM THE FIRST SUBSTRUCTURE IN THE ARRAY OLDNMS
!     - NEWNM RESULTS FROM COMBINING THE FIRST I SUBSTRUCTURES IN THE
!       ARRAY OLDNMS IF NUMB = I
 
!     THE OUTPUT VARIABLE ITEST TAKES ON ONE OF THE FOLLOWING VALUES
!          4  IF ONE  OR MORE SUBSTRUCTURES IN OLDNMS DO NOT EXIST
!          7  IF NEWNM ALREADY EXISTS
!          8  IF ONE OF THE SUBSTRUCTURES IN OLDNMS HAS ALREADY
!             BEEN USED IN A REDUCTION OR COMBINATION
!          1  OTHERWISE
 
!     IF ITEST IS SET TO 4, NUMB WILL BE SET TO THE NUMBER OF
!     SUBSTRUCTURES IN OLDNMS THAT DO NOT EXIST AND THE FIRST NUMB NAMES
!     IN OLDNMS WILL BE SET TO THE NAMES OF THOSE SUBSTRUCTURES THAT DO
!     NOT EXIST.  BIT IBIT OF THE FIRST MDI WORD IS SET TO INDICATE THE
!     APPROPRIATE TYPE OF SUBSTRUCTURE. IF IBIT IS ZERO NO CHANGE IS
!     MADE TO THE MDI
 
 
 INTEGER, INTENT(IN OUT)                  :: newnm(2)
 INTEGER, INTENT(IN OUT)                  :: numb
 INTEGER, INTENT(OUT)                     :: oldnms(14)
 INTEGER, INTENT(OUT)                     :: itest
 INTEGER, INTENT(IN OUT)                  :: ibit
 EXTERNAL        lshift,andf,orf,complf
 LOGICAL :: ditup,mdiup
 DIMENSION  iold(7),nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     iodum(8),mdi,mdipbn,mdilbn,mdibl, nxtdum(15),ditup,mdiup
 DATA    iempty/ 4H    /, nmsbr / 4HSETL,4HVL  /
 DATA    ll,cs , hl    /  2,2,2 /
 DATA    ib    / 1     /
 
 CALL chkopn (nmsbr(1))
 itest = 1
 CALL fdsub (newnm(1),i)
 IF (i  /=  -1) GO TO 500
 IF (numb == 0) GO TO 20
 
!     MAKE SURE THAT ALL THE SUBSTRUCTURES IN OLDNMS DO EXIST.
 
 icount = 0
 
 DO  i = 1,numb
   k = 2*(i-1) + 1
   CALL fdsub (oldnms(k),iold(i))
   IF (iold(i) > 0) CYCLE
   icount = icount + 1
   kk = 2*(icount-1) + 1
   oldnms(kk  ) = oldnms(k  )
   oldnms(kk+1) = oldnms(k+1)
 END DO
 
 IF (icount == 0) GO TO 20
 numb = icount
 GO TO 510
 20 CALL crsub (newnm(1),inew)
 IF (numb == 0) RETURN
 
!     NEWNM IS NOT A BASIC SUBSTRUCTURE (LEVEL 0).
!     UPDATE NEWNM S DIRECTORY IN THE MDI.
 
 CALL fmdi (inew,imdi)
 llmask = complf(lshift(1023,20))
 buf(imdi+ll) = orf(andf(buf(imdi+ll),llmask),lshift(iold(1),20))
 IF (ibit /= 0) buf(imdi+ib) = orf(buf(imdi+ib),lshift(1,ibit))
 mdiup = .true.
 
!     UPDATE IN THE MDI THE DIRECTORIES OF THE SUBSTRUCTURES IN OLDNMS.
 
 IF (numb > 7) numb = 7
 maskcs = complf(lshift(1023,10))
 
 DO  i = 1,numb
   CALL fmdi (iold(i),imdi)
   IF (andf(buf(imdi+hl),1023) == 0) GO TO 40
   icount = i
   GO TO 520
   40 buf(imdi+hl) = orf(buf(imdi+hl),inew)
   mdiup = .true.
   IF (numb == 1) RETURN
   IF (i == numb) EXIT
   buf(imdi+cs) = orf(andf(buf(imdi+cs),maskcs),lshift(iold(i+1),10))
 END DO
 
 130 buf(imdi+cs) = orf(andf(buf(imdi+cs),maskcs),lshift(iold(1),10))
 RETURN
 
!     NEWNM ALREADY EXISTS.
 
 500 itest = 7
 RETURN
 
!     ONE OR MORE OF THE SUBSTRUCTURES IN OLDNMS DO NOT EXIST.
 
 510 itest = 4
 RETURN
 
!     ONE OF THE SUBSTRUCTURES IN OLDNMS HAS ALREADY BEEN USED IN A
!     REDUCTION OR COMBINATION.  REMOVE ALL CHANGES THAT HAVE BEEN MADE.
 
 520 itest = 8
 CALL fdit (inew,idit)
 buf(idit  ) = iempty
 buf(idit+1) = iempty
 ditup = .true.
 IF (2*inew /= ditsiz) GO TO 525
 ditsiz = ditsiz - 2
 525 ditnsb = ditnsb - 1
 CALL fmdi (inew,imdi)
 buf(imdi+ll) = andf(buf(imdi+ll),llmask)
 mdiup  = .true.
 icount = icount - 1
 IF (icount < 1) RETURN
 
 DO  i = 1,icount
   CALL fmdi (iold(i),imdi)
   buf(imdi+hl) = andf(buf(imdi+hl),complf(1023))
   buf(imdi+cs) = andf(buf(imdi+cs),maskcs)
   mdiup = .true.
 END DO
 
 RETURN
END SUBROUTINE setlvl
