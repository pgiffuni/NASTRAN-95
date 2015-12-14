SUBROUTINE dstroy (NAME,itest,image,imore,lim)
     
!     DESTROYS THE SUBSTRUCTURE NAME BY DELETING ITS DIRECTORY FROM THE
!     MDI AND ITS NAME FROM THE DIT.  NO OPERATION WILL TAKE PLACE IF
!     NAME IS AN IMAGE SUBSTRUCTURE.  IF NAME IS A SECONDARY SUBSTRUC-
!     TURE, IT IS DELETED FROM THE LIST OF SECONDARY SUBSTRUCTURES TO
!     WHICH IT BELONGS, AND ITS IMAGE CONTRIBUTING TREE IS DESTROYED.
!     IF NAME IS A PRIMARY SUBSTRUCTURE, ALL ITS SECONDARY SUBSTRUCTURES
!     ARE ALSO DESTROYED.  IN ALL CASES, ALL THE SUBSTRUCTURES DERIVED
!     FROM THE SUBSTRUCTURE BEING DESTROYED ARE ALSO DESTROYED, AND
!     CONNECTIONS WITH OTHER SUBSTRUCTURES ARE DELETED.
 
!     THE BLOCKS OCCUPIED BY THE ITEM ARE RETURNED TO THE LIST OF FREE
!     BLOCKS IF THEY BELONG TO THE SPECIFIED SUBSTRUCTURE
 
!     THE OUTPUT VARIABLE ITEST TAKES ONE OF THE FOLLOWING VALUES.
!        1  NORMAL RETURN
!        4  IF NAME DOES NOT EXIST
!        6  IF NAME IS AN IMAGE SUBSTRUCTURE
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(OUT)                     :: itest
 INTEGER, INTENT(OUT)                     :: image(1)
 INTEGER, INTENT(OUT)                     :: imore(1)
 INTEGER, INTENT(IN OUT)                  :: lim
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: ditup,mdiup
 INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     mdi,mdipbn,mdilbn,mdibl,blksiz,dirsiz,ps,ss,is,  &
     ll,cs,hl,andf,orf,rshift,complf
 DIMENSION  nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     iodum(8),mdi,mdipbn,mdilbn,mdibl, nxtdum(15),ditup,mdiup
 COMMON /sys   / blksiz,dirsiz,sys(3),ifrst
 COMMON /itemdt/ nitem,item(7,1)
 DATA    ps,ss,  is,ll,cs,hl / 1,1,1,2,2,2    /
 DATA    iempty/ 4H    /
 DATA    indsbr/ 3     /, nmsbr /4HDSTR,4HOY  /
 
 CALL chkopn (nmsbr(1))
 itest = 1
 itop  = 0
 imtop = 0
 CALL fdsub (NAME(1),INDEX)
 IF (INDEX == -1) GO TO 1000
 maskm = complf(lshift(1023,10))
 maskl = complf(lshift(1023,20))
!                           1023 = 2**10 - 1
 
!     SAVE ALL CONNECTIONS WITH OTHER SUBSTRUCTURES.
 
 10 CALL fmdi (INDEX,imdi)
 20 i     = buf(imdi+ps)
 indps = andf(i,1023)
 indss = rshift(andf(i,1048575),10)
!                           1048575 = 2**20 - 1
 indis = andf(i,1073741824)
!                    1073741824 = 2**30
 i     = buf(imdi+ll)
 indhl = andf(i,1023)
 indcs = rshift(andf(i,1048575),10)
 indll = rshift(andf(i,1073741823),20)
!                           1073741823 = 2**30 - 1
 IF (indis > 0) GO TO 1010
 IF (indps == 0) GO TO 60
 ASSIGN 30 TO iret1
 GO TO 300
 
!     REMOVE INDEX FROM THE LIST OF SUBSTRUCTURES THAT ARE SECONDARY TO
!     INDPS.
 
 30 isave = indps
 40 CALL fmdi (isave,imdi)
 isave = rshift(andf(buf(imdi+ss),1048575),10)
 IF (isave == 0) GO TO 50
 IF (isave /= INDEX) GO TO 40
 buf(imdi+ss) = orf(andf(buf(imdi+ss),maskm),lshift(indss,10))
 mdiup = .true.
 IF (indll == 0) GO TO 120
 ill   = indll
 indll = 0
 isave = INDEX
 50 ASSIGN 120 TO iret2
 GO TO 330
 
!     PRIMARY SUBSTRUCTURE.
!     RETURN THE BLOCKS USED BY ALL ITEMS TO THE LIST OF FREE BLOCKS.
 
 60 DO  j = ifrst,dirsiz
   ibl = andf(buf(imdi+j),65535)
!                            65535 = 2**16 - 1
   IF (ibl > 0 .AND. ibl /= 65535) CALL retblk (ibl)
 END DO
 IF (indss == 0) GO TO 130
 
!     THE PRIMARY SUBSTRUCTURE BEING DESTROYED HAS SECONDARY EQUIVALENT
!     SUBSTRUCTURES.  MUST DESTROY ALL OF THEM.
 
 ASSIGN 320 TO iret1
 ASSIGN 90  TO iret2
 isv   = indss
 80 isave = isv
 CALL fmdi (isave,imdi)
 isv = rshift(andf(buf(imdi+ss),1048575),10)
 iis = andf(buf(imdi+is),1073741824)
 IF (iis > 0) GO TO 110
 
!     THE SECONDARY SUBSTRUCTURE IS NOT AN IMAGE SUBSTRUCTURE.  ADD ITS
!     INDEX TO THE LIST (IMORE) OF SUBSTRUCTURES TO BE DESTROYED LATER.
 
 itop = itop + 1
 IF (itop > lim) GO TO 1030
 imore(itop) = isave
 GO TO 300
 
!     UPDATE THE MDI OF THE SECONDARY SUBSTRUCTURE WITH INDEX ISAVE.
 
 90 CALL fmdi (isave,imdi)
 buf(imdi+ps) = 0
 buf(imdi+ll) = andf(buf(imdi+ll),maskl)
 DO  j = ifrst,dirsiz
   buf(imdi+j) = 0
 END DO
 mdiup = .true.
 110 IF (isv /= 0) GO TO 80
 
!     BACK TO THE SUBSTRUCTURE WITH INDEX  INDEX .
!     DELETE ITS DIRECTORY FROM THE MDI.
 
 120 CALL fmdi (INDEX,imdi)
 130 DO  j = 1,dirsiz
   buf(imdi+j) = 0
 END DO
 mdiup = .true.
 
!     DELETE SUBSTRUCTURE NAME FROM THE DIT.
 
 CALL fdit (INDEX,jdit)
 buf(jdit  ) = iempty
 buf(jdit+1) = iempty
 ditup = .true.
 IF (INDEX*2 /= ditsiz) GO TO 150
 ditsiz = ditsiz - 2
 150 ditnsb = ditnsb - 1
 IF (indcs == 0) GO TO 180
 
!     DELETE LINK THROUGH COMBINED SUBSTRUCTURES, AND REMOVE ITEMS
!     CREATED AS A RESULTS OF THE COMBINE OR REDUCE.
!     THESE ITEMS WILL BE RETURNED TO THE LIST OF FREE BLOCKS.
 
 160 IF (indcs == INDEX) GO TO 180
 CALL fmdi (indcs,imdi)
 indcs = rshift(andf(buf(imdi+cs),1048575),10)
 173 buf(imdi+hl) = andf(buf(imdi+hl),complf(1023))
 buf(imdi+cs) = andf(buf(imdi+cs),maskm)
 DO  j = 1,nitem
   IF (item(6,j) == 0) CYCLE
   itm = j + ifrst - 1
   ibl = andf(buf(imdi+itm),65535)
   IF (ibl > 0 .AND. ibl /= 65535) CALL retblk (ibl)
   buf(imdi+itm) = 0
 END DO
 mdiup = .true.
 IF (indcs == 0) GO TO 1020
 GO TO 160
 180 IF (indll == 0) GO TO 190
 
!     SUBSTRUCTURE WAS THE RESULT OF COMBINING LOWER LEVEL SUBSTRUCTURES
!     TOGETHER.  UPDATE THE MDI ACCORDINGLY.
 
 CALL fmdi (indll,imdi)
 indcs = rshift(andf(buf(imdi+cs),1048575),10)
 INDEX = indll
 indll = 0
 IF (indcs == 0) indcs = INDEX
 GO TO 173
 190 IF (indhl == 0) GO TO 220
 
!     A HIGHER LEVEL SUBSTRUCTURE WAS DERIVED FROM THE ONE BEING
!     DESTROYED. DESTROY THE HIGHER LEVEL SUBSTRUCTURE.
 
 INDEX = indhl
 CALL fmdi (INDEX,imdi)
 buf(imdi+ll) = andf(buf(imdi+ll),maskl)
 mdiup = .true.
 GO TO 20
 220 IF (itop == 0) RETURN
 
!     MORE SUBSTRUCTURES TO DESTROY.
 
 INDEX = imore(itop)
 itop  = itop - 1
 GO TO 10
 
!     INTERNAL SUBROUTINE.
!     RETURN TO THE LIST OF FREE BLOCKS THE BLOCKS USED BY A
!     SECONDARY SUBSTRUCTURE.
!     THESE BLOCKS INCLUDE THE FOLLOWING ITEMS
 
!     ITEMS COPIED DURING A EQUIV OPERATION
!     SOLUTION ITEMS
!     ITEMS PRODUCED BY A COMBINE OR REDUCE OPERATION
 
 300 DO  j = 1,nitem
   IF (item(5,j) == 0) CYCLE
   itm = j + ifrst - 1
   ibl = andf(buf(imdi+itm),65535)
   IF (ibl > 0 .AND. ibl /= 65535) CALL retblk (ibl)
   buf(imdi+itm) = 0
 END DO
 GO TO iret1, (30,320)
 
!     INTERNAL SUBROUTINE.
!     BUILD A LIST IMAGE OF ALL THE IMAGE SUBSTRUCTURES CONTRIBUTING TO
!     THE SECONDARY SUBSTRUCTURE WITH INDEX ISAVE, AND DELETE EACH IMAGE
!     SUBSTRUCTURE FROM THE LIST OF SECONDARY SUBSTRUCTURES TO WHICH IT
!     BELONGS.
 
 320 CALL fmdi (isave,imdi)
 ill = rshift(andf(buf(imdi+ll),1073741823),20)
 IF (ill == 0) GO TO iret2, (90,120)
 330 imtop = 1
 image(imtop) = ill
 icount = 1
 ihere  = image(icount)
 350 CALL fmdi (ihere,imdi)
 i   = buf(imdi+ps)
 ips = andf(i,1023)
 iss = rshift(andf(i,1048575),10)
 iis = andf(i,1073741824)
 i   = buf(imdi+ll)
 ill = rshift(andf(i,1073741823),20)
 ics = rshift(andf(i,1048575),10)
 IF (iis == 0) GO TO 1010
 
!     DELETE THE SUBSTRUCTURE WITH INDEX IHERE FROM THE MDI AND THE DIT.
!     RETURN THE BLOCKS USED BY THE IMAGE SUBSTRUCTURE TO THE LIST OF
!     FREE BLOCKS.  THIS INCLUDES THE FOLLOWING ITEMS
 
!     ITEMS COPIED DURING A EQUIV OPERATION
!     SOLUTION ITEMS
 
 DO  j = 1,nitem
   IF (item(4,j) == 0) CYCLE
   itm = j + ifrst - 1
   ibl = andf(buf(imdi+itm),65535)
   IF (ibl > 0 .AND. ibl /= 65535) CALL retblk (ibl)
   buf(imdi+itm) = 0
 END DO
 DO  j = 1,dirsiz
   buf(imdi+j) = 0
 END DO
 mdiup = .true.
 CALL fdit (ihere,idit)
 buf(idit  ) = iempty
 buf(idit+1) = iempty
 ditup = .true.
 IF (ihere*2 /= ditsiz) GO TO 370
 ditsiz = ditsiz - 2
 370 ditnsb = ditnsb - 1
 
!     DELETE POINTERS TO IHERE.
 
 icheck = ips
 380 CALL fmdi (icheck,imdi)
 icheck = rshift(andf(buf(imdi+ss),1048575),10)
 IF (icheck == 0) GO TO 390
 IF (icheck /= ihere) GO TO 380
 buf(imdi+ss) = orf(andf(buf(imdi+ss),maskm),lshift(iss,10))
 mdiup = .true.
 
!     ARE THERE MORE SUBSTRUCTURES TO ADD TO THE LIST IMAGE
 
 390 IF (ill == 0) GO TO 410
 DO  j = 1,imtop
   IF (image(j) == ill) GO TO 410
 END DO
 imtop = imtop + 1
 image(imtop) = ill
 410 IF (ics == 0) GO TO 430
 DO  j = 1,imtop
   IF (image(j) == ics) GO TO 430
 END DO
 imtop = imtop + 1
 IF (imtop > lim) GO TO 1030
 image(imtop) = ics
 
!     ARE THERE MORE SUBSTRUCTURES ON THE LIST IMAGE
 
 430 IF (icount == imtop) GO TO iret2, (90,120)
 icount = icount + 1
 ihere  = image(icount)
 GO TO 350
 
!     NAME DOES NOT EXIST.
 
 1000 itest = 4
 RETURN
 
!     NAME IS AN IMAGE SUBSTRUCTURE.
 
 1010 itest = 6
 RETURN
 
!     ERROR MESSAGES.
 
 1020 CALL errmkn (indsbr,8)
 1030 CALL mesage (-8,0,nmsbr)
 RETURN
END SUBROUTINE dstroy
