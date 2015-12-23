SUBROUTINE seteq (name1,name2,prefx,dry2,itest,imore,lim)
     
!     SETS THE SUBSTRUCTURE NAME2 EQUIVALENT TO THE SUBSTRUCTURE NAME1.
!     THE OUTPUT VARIABLE ITEST TAKES ON ONE OF THE FOLLOWING VALUES
 
!         4  IF NAME1 DOES NOT EXIST
!         8  IF DRY DOES NOT EQUAL ZERO AND NAME2 OR ONE OF THE NEW
!            NAMES ALREADY EXISTS
!         9  IF DRY IS EQUAL TO ZERO AND NAME2 OR ONE OF THE NEW NAMES
!            DOES NOT EXIST
!         1  OTHERWISE
 
 IMPLICIT INTEGER (a-z) 
 INTEGER, INTENT(IN)                      :: name1(2)
 INTEGER, INTENT(IN)                      :: name2(2)
 INTEGER, INTENT(IN OUT)                  :: prefx
 INTEGER, INTENT(IN)                      :: dry2
 INTEGER, INTENT(OUT)                     :: itest
 INTEGER, INTENT(OUT)                     :: imore(1)
 INTEGER, INTENT(IN OUT)                  :: lim
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: ditup,mdiup,more
 DIMENSION  isave(50),namnew(2),  nmsbr(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     io ,iodum(7),mdi,mdipbn,mdilbn,mdibl, nxtdum(15),ditup,mdiup
 COMMON /sys   / blksiz,dirsiz,sys(3),ifrst
 COMMON /output/ title(96),subtit(96)
 COMMON /system/ nbuff,nout,dum(36),nbpc,nbpw,ncpw
 COMMON /itemdt/ nitem,item(7,1)
 DATA    ps, ss, ib, ll, cs, hl, bb,   ird, iwrt,  indsbr  /  &
     1 , 1 , 1 ,  2,  2,  2,  1,     1,    2,      15  /
 DATA    iempty, mask,   nmsbr         / 4H    , 4HMASK, 4HSETE,4HQ    /
 
 CALL chkopn (nmsbr(1))
 IF (nitem+ifrst-1 > 50) GO TO 970
 dry   = dry2
 itest = 1
 CALL fdsub (name1(1),ind1)
 IF (ind1 == -1) GO TO 900
 mask   = andf(mask,2**(nbpw-4*nbpc)-1)
 maskss = complf(lshift(1023,10))
 maskll = complf(lshift(1023,20))
 maskbb = lshift(1023,20)
 
!     IF NAME2 EXISTS - VERIFY THAT IT IS MARKED EQUIVALENT TO NAME1.
!     NAME2 MAY ALREADY EXIST FOR RUN=GO OR OPTIONS=PA
 
 CALL fdsub (name2(1),ind2)
 IF (ind2 == -1) GO TO 10
 dry = 0
 
 CALL fmdi (ind2,imdi)
 ips = andf(buf(imdi+ps),1023)
 IF (ips ==    0) GO TO 920
 IF (ips == ind1) GO TO 10
 CALL fmdi (ind1,imdi)
 ipp = andf(buf(imdi+ps),1023)
 IF (ips /= ipp) GO TO 920
 
!     STEP 1.  MAKE A LIST OF ALL THE SUBSTRUCTURES CONTRIBUTING TO THE
!     SUBSTRUCTURE NAME1, AND STORE IT IN THE ARRAY IMORE
 
 10 itop  = 1
 imore(itop) = ind1
 iptr  = 1
 20 CALL fmdi (ind1,imdi)
 i     = buf(imdi+ll)
 indll = rshift(andf(i,1073741823),20)
 indcs = rshift(andf(i,1048575)   ,10)
 IF (indll == 0) GO TO 40
 DO  j = 1,itop
   IF (imore(j) == indll) GO TO 40
 END DO
 itop  = itop + 1
 IF (itop > lim) GO TO 960
 imore(itop) = indll
 40 IF (indcs == 0 .OR. iptr == 1) GO TO 60
 DO  j = 1,itop
   IF (imore(j) == indcs) GO TO 60
 END DO
 itop  = itop + 1
 IF (itop > lim) GO TO 960
 imore(itop) = indcs
 60 IF (iptr == itop) GO TO 100
 iptr  = iptr + 1
 ind1  = imore(iptr)
 GO TO 20
 
!     STEP 2.  CREATE AN IMAGE SUBSTRUCTURE FOR EACH SUBSTRUCTURE IN THE
!     ARRAY IMORE, AND STORE ITS INDEX IN THE ARRAY IMAGE.  NOTE THAT
!     SINCE IMORE(1) CONTAINS THE INDEX OF NAME1, IMAGE(1) WILL CONTAIN
!     THE INDEX OF NAME2
!     FOR EACH NEW NAME CHECK THAT MAKING ROOM FOR THE PREFIX DOES NOT
!     TRUNCATE THE NAME
 
 100 IF (iptr /= 1) GO TO 110
 CALL fdsub (name2(1),i)
 GO TO 120
 110 CALL fdit (ind1,idit)
 first = klshft(krshft(prefx,ncpw-1),ncpw-1)
 rest  = klshft(krshft(buf(idit),ncpw-3),ncpw-4)
 namnew(1) = orf(orf(first,rest),mask)
 first = klshft(krshft(buf(idit),ncpw-4),ncpw-1)
 rest  = klshft(krshft(buf(idit+1),ncpw-3),ncpw-4)
 namnew(2)= orf(orf(first,rest),mask)
 IF (khrfn1(iempty,4,buf(idit+1),4) /= iempty)  &
     WRITE (nout,850) uwm,namnew,buf(idit),buf(idit+1)
 CALL fdsub (namnew(1),i)
 120 IF (dry /= 0) GO TO 130
 IF (i   /= -1) GO TO 170
 GO TO 910
 130 IF (i == -1) GO TO 150
 iptr = iptr + 1
 IF (iptr > itop) GO TO 920
 DO  i = iptr,itop
   image = imore(lim+i)
   CALL fdit (image,idit)
   buf(idit  ) = iempty
   buf(idit+1) = iempty
   ditup = .true.
 END DO
 GO TO 920
 150 IF (iptr /= 1) GO TO 160
 CALL crsub (name2(1),i)
 GO TO 170
 160 CALL crsub (namnew(1),i)
 170 imore(iptr+lim) = i
 IF (iptr == 1) GO TO 200
 iptr = iptr - 1
 ind1 = imore(iptr)
 GO TO 100
 
!     STEP 3.  BUILD THE MDI OF NAME2, AND OF ALL IMAGE SUBSTRUCTURES
 
 200 ind2 = i
 210 CALL fmdi (ind1,imdi)
 DO  j = 1,dirsiz
   isave(j) = buf(imdi+j)
 END DO
 
!     SET THE SS ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND1
 
 IF (dry == 0) GO TO 230
 buf(imdi+ss) = orf(andf(buf(imdi+ss),maskss),lshift(ind2,10))
 mdiup = .true.
 230 CALL fmdi (ind2,imdi)
 IF (dry == 0) GO TO 420
 i = isave(ps)
 
!     SET THE PS ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 ips = andf(i,1023)
 IF (ips == 0) GO TO 240
 buf(imdi+ps) = ips
 GO TO 250
 240 buf(imdi+ps) = ind1
 
!     SET THE SS ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 250 iss = rshift(andf(i,1048575),10)
 IF (iss == 0) GO TO 260
 buf(imdi+ss) = orf(andf(buf(imdi+ss),maskss),lshift(iss,10))
 
!     SET THE BB ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 260 ibs = andf(i,maskbb)
 buf(imdi+bb) = orf(andf(buf(imdi+bb),maskll),ibs)
 i = isave(ll)
 
!     SET THE HL ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 IF (iptr == 1) GO TO 300
 ihl = andf(i,1023)
 IF (ihl == 0) GO TO 280
 ASSIGN 270 TO iret
 iwant = ihl
 GO TO 320
 270 buf(imdi+hl) = ifnd
 
!     SET THE CS ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 280 ics = rshift(andf(i,1048575),10)
 IF (ics == 0) GO TO 300
 ASSIGN 290 TO iret
 iwant = ics
 GO TO 320
 290 buf(imdi+cs) = orf(andf(buf(imdi+cs),maskss),lshift(ifnd,10))
 
!     SET THE LL ENTRY FOR THE SUBSTRUCTURE WITH INDEX IND2
 
 300 ill = rshift(andf(i,1073741823),20)
 IF (ill == 0) GO TO 400
 ASSIGN 310 TO iret
 iwant = ill
 GO TO 320
 310 buf(imdi+ll) = orf(andf(buf(imdi+ll),maskll),lshift(ifnd,20))
 GO TO 400
 
!     FIND THE INDEX OF THE IMAGE SUBSTRUCTURE TO THE SUBSTRUCTURE WITH
!     INDEX IWANT.  STORE THE FOUND INDEX IN IFND
 
 320 DO  k = 1,itop
   IF (imore(k) /= iwant) CYCLE
   ifnd = imore(lim+k)
   GO TO iret, (270,290,310)
 END DO
 GO TO 930
 
!     SET THE POINTERS OF THE ITEMS BELONGING TO THE SUBSTRUCTURE WITH
!     INDEX IND2
 
 400 DO  j = ifrst,dirsiz
   buf(imdi+j) = 0
 END DO
 420 IF (iptr == 1) GO TO 440
 
!     IMAGE SUBSTRUCTURE - SET POINTERS TO SHARED ITEMS AND SET IB BIT
 
 DO  j = 1,nitem
   IF (item(4,j) /= 0) CYCLE
   itm = j + ifrst - 1
   IF (buf(imdi+itm) == 0) buf(imdi+itm) = isave(itm)
 END DO
 buf(imdi+ib) = orf(buf(imdi+ib),lshift(1,30))
 GO TO 500
 
!     SECONDARY SUBSTRUCTURE - SET POINTERS TO SHARED ITEMS
 
 440 DO  j = 1,nitem
   IF (item(5,j) /= 0) CYCLE
   itm = j + ifrst - 1
   IF (buf(imdi+itm) == 0) buf(imdi+itm) = isave(itm)
 END DO
 
!     COPY APPROPRIATE ITEMS OF NAME1 AND WRITE THEM FOR
!     NAME2 AFTER CHANGING NAME1 TO NAME2 AND INSERTING THE NEW PREFIX
!     TO THE NAMES OF ALL CONTRIBUTING SUBSTRUCTURES
 
 500 DO  j = 1,nitem
   IF (item(3,j) == 0) CYCLE
   kk = j + ifrst - 1
   IF (buf(imdi+kk) /= 0) CYCLE
   irdbl = andf(isave(kk),jhalf)
   IF (irdbl /= 0 .AND. irdbl /= jhalf) GO TO 510
   buf(imdi+kk) = isave(kk)
   CYCLE
   510 CALL sofio (ird,irdbl,buf(io-2))
   CALL fdit (ind2,idit)
   buf(io+1) = buf(idit  )
   buf(io+2) = buf(idit+1)
   CALL getblk (0,iwrtbl)
   IF (iwrtbl == -1) GO TO 940
   newblk = iwrtbl
   numb = item(3,j)/1000000
   MIN  = (item(3,j) - numb*1000000)/1000
   inc  = item(3,j) - numb*1000000 - MIN*1000
   numb = buf(io+numb)
   IF (numb > 1 .OR. ill /= 0 .OR. iptr /= 1) GO TO 530
   
!     BASIC SUBSTRUCTURE
   
   buf(io+MIN  ) = name2(1)
   buf(io+MIN+1) = name2(2)
   more = .false.
   GO TO 580
   
!     NOT A BASIC SUBSTRUCTURE
   
   530 IF (numb <= (blksiz-MIN+1)/inc) GO TO 540
   numb = numb - (blksiz-MIN+1)/inc
   MAX  = blksiz
   more = .true.
   GO TO 550
   540 MAX  = MIN + inc*numb - 1
   more = .false.
   
!     INSERT THE NEW PREFIX TO THE NAMES OF ALL CONTRIBUTING SUBSTRUC-
!     TURES
!     IF THE COMPONENT IS FOR MODAL DOF ON THE SECONDARY SUBSTRUCTURE,
!     USE THE ACTUAL NAME INSTEAD OF ADDING A PREFIX
   
   550 DO  k = MIN,MAX,inc
     IF (buf(io+k) == name1(1) .AND. buf(io+k+1) == name1(2)) GO TO 560
     first = klshft(krshft(prefx,ncpw-1),ncpw-1)
     rest  = klshft(krshft(buf(io+k  ),ncpw-3),ncpw-4)
     first2= klshft(krshft(buf(io+k  ),ncpw-4),ncpw-1)
     rest2 = klshft(krshft(buf(io+k+1),ncpw-3),ncpw-4)
     buf(io+k  ) = orf(orf(first ,rest ),mask)
     buf(io+k+1) = orf(orf(first2,rest2),mask)
     CYCLE
     
     560 buf(io+k  ) = name2(1)
     buf(io+k+1) = name2(2)
   END DO
   
!     WRITE OUT UPDATED DATA BLOCK
   
   580 CALL sofio (iwrt,iwrtbl,buf(io-2))
   CALL fnxt (irdbl,inxt)
   IF (MOD(irdbl,2) == 1) GO TO 590
   next = andf(rshift(buf(inxt),ihalf),jhalf)
   GO TO 600
   590 next = andf(buf(inxt),jhalf)
   600 IF (next == 0) GO TO 620
   
!     MORE BLOCKS TO COPY
   
   irdbl = next
   CALL getblk (iwrtbl,next)
   IF (next /= -1) GO TO 610
   CALL retblk (newblk)
   GO TO 940
   610 iwrtbl = next
   CALL sofio (ird,irdbl,buf(io-2))
   MIN = 1
   IF (more) GO TO 530
   GO TO 580
   
!     NO MORE BLOCKS TO COPY.  UPDATE MDI OF NAME2
   
   620 buf(imdi+kk) = orf(lshift(rshift(isave(kk),ihalf),ihalf),newblk)
 END DO
 
 mdiup = .true.
 IF (iptr == itop) GO TO 720
 iptr = iptr + 1
 ind1 = imore(iptr    )
 ind2 = imore(iptr+lim)
 GO TO 210
 
!     WRITE USER MESSAGES
 
 720 IF(dry == 0) GO TO 780
 DO  i = 1,96
   subtit(i) = iempty
 END DO
 CALL page
 CALL page2 (-4)
 WRITE (nout,800) name2,name1
 image = imore(lim+1)
 CALL fmdi (image,imdi)
 ips = andf(buf(imdi+1),1023)
 CALL fdit (ips,i)
 CALL page2 (-2)
 WRITE (nout,810) name2,buf(i),buf(i+1)
 iptr = 2
 IF (iptr > itop) GO TO 990
 CALL page2 (-2)
 WRITE (nout,820)
 740 DO  i = 1,16
   imore(i) = iempty
 END DO
 j = 1
 760 image = imore(iptr+lim)
 CALL fdit (image,i)
 imore(j  ) = buf(i  )
 imore(j+1) = buf(i+1)
 iptr = iptr + 1
 IF (iptr > itop) GO TO 770
 j = j + 2
 IF (j < 16) GO TO 760
 770 CALL page2 (-2)
 WRITE (nout,830) (imore(j),j=1,16)
 IF (iptr <= itop) GO TO 740
 GO TO 990
 
!     DRY RUN - PRINT MESSAGE INDICATING ONLY ADDITIONS MADE
 
 780 CALL page2 (-3)
 WRITE (nout,840) uim,name2,name1,name2
 GO TO 990
 
 800 FORMAT (32X,67HS u b s t r u c t u r e   e q u i v a l e n c e   o &
     &p e r a t i o n ,///23X,13HSUBSTRUCTURE ,2A4,56H has been created &
     &and marked equivalent TO substructure ,2A4)
 810 FORMAT (1H0,22X,28HTHE primary substructure of ,2A4,4H is ,2A4)
 820 FORMAT (1H0,22X, 56HTHE following image substructures have been ge&
     &nerated --)
 830 FORMAT (1H0,22X,10(2A4,2X))
 840 FORMAT (a29,' 6228, SUBSTRUCTURE ',2A4,' IS ALREADY AN EQUIVALENT'  &
     ,      ' SUBSTRUCTURE TO ',2A4, /36X,'ONLY ITEMS NOT PREVIOUSLY ',  &
     'EXISTING FOR ',2A4,' HAVE BEEN MADE EQUIVALENT.')
 850 FORMAT (a25,' 6236, DURING THE CREATION OF A NEW IMAGE SUBSTRUC',  &
     'TURE NAMED ',2A4,' THE LAST CHARACTER ', /5X,  &
     'OF SUBSTRUCTURE NAMED ',2A4,' WAS TRUNCATED TO MAKE ROOM',  &
     ' FOR THE PREFIX.')
 
!     ERROR CONDITIONS
 
 900 itest = 4
 GO TO 990
 910 itest = 9
 GO TO 990
 920 itest = 8
 GO TO 990
 930 CALL errmkn (indsbr,3)
 940 WRITE  (nout,950) ufm
 950 FORMAT (a23,' 6223, SUBROUTINE SETEQ - THERE ARE NO MORE FREE ',  &
     'BLOCKS AVAILABLE ON THE SOF.')
 k = -37
 GO TO 980
 960 k = -8
 GO TO 980
 970 CALL errmkn (indsbr,10)
 980 CALL sofcls
 CALL mesage (k,0,nmsbr)
 
 990 RETURN
END SUBROUTINE seteq
