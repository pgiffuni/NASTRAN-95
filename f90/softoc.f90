SUBROUTINE softoc
     
!     SOF TABLE OF CONTENTS ROUTINE
 
 
!     THE CURRENT SUBSTRUCTURE TYPE BIT POSITIONS ARE -
 
!        NO BIT - BASIC SUBSTRUCTURE (EXCEPT IMAGE BIT)
!        BIT 30 - IMAGE SUBSTRUCTURE
!            29 - COMBINED SUBSTRUCTURE
!            28 - GUYAN REDUCTUION SUBSTRUCTURE
!            27 - MODAL REDUCTION SUBSTRUCTURE
!            26 - COMPLEX MODAL REDUCTION SUBSTRUCTURE
 
!     TO ADD A NEW SUBSTRUCTURE TYPE BIT THE FOLLOWING UPDATES ARE
!     REQUIRED.
 
!        1) INCREASE THE DEMENSION OF TYPE.
!        2) INCREASE THE VALUE OF NTYPE IN THE DATA STATEMENT.
!        3) ADD A NEW BCD TYPE VALUE TO THE DATA STATEMENT.
 
 
!     THIS ROUTINE IS CURRENTLY CODED TO HANDLE UP TO 27 SOF ITEMS
!     AUTOMATICALLY.
!     TO INCREASE THIS TO 40 ITEMS PERFORM THE FOLLOWING UPDATES.
 
!        1) CHANGE THE DIMENSION OF HDR TO (40,4)
!        2) CHANGE THE DIMENSION OF ITM TO (40)
!        3) CHANGE THE VALUE OF MAXITM IN THE DATA STATEMENT TO 40
!        4) CHANGE THE INNER GROUPS ON FORMAT 80 TO 39(A1,1X),A1
!        5) CHANGE THE INNER GROUP ON FORMAT 100 TO 39(A1,1X),A1
 
 EXTERNAL        lshift,rshift,andf
 INTEGER :: avblks,BLANK,ditnsb,buf,ssname(2),andf,ss,ps,cs,  &
     hl,rshift,dirsiz,sofsiz,ditsiz,num(10),blksiz,  &
     hiblk,filsiz,TYPE(5),itm(27),hdr(27,4)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /machin/ mach,ihalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl
 COMMON /sys   / blksiz,dirsiz,supsiz,avblks,hiblk,ifrst
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10)
 COMMON /system/ sysbuf,nout,z1(6),nlpp,y(2),line,z2(26),nbpc,nbpw
 COMMON /itemdt/ nitem,item(7,1)
 DATA    TYPE  / 2HB , 2HC , 2HR , 2HM , 2HCM /
 DATA    num   / 1H1, 1H2, 1H3, 1H4, 1H5, 1H6 ,1H7, 1H8, 1H9, 1H0 /
 DATA    BLANK / 4H     /
 DATA    image / 4HI    /
 DATA    ntype / 6      /
 DATA    maxitm/ 27     /
 
 nitm = nitem
 IF (nitm <= maxitm) GO TO 10
 nitm = maxitm
 WRITE  (nout,6237) swm,maxitm
 6237 FORMAT (a27,' 6237, THE SOFTOC ROUTINE CAN HANDLE ONLY',i4,  &
     ' ITEMS.', /34X,'ADDITIONAL ITEMS WILL NOT BE SHOWN')
 
!     SET UP HEADINGS AND MASKS
 
 10 nshft = 0
 DO  i = 1,4
   DO  j = 1,nitm
     hdr(j,i) = klshft(item(1,j),nshft/nbpc)
   END DO
   k = nitm + 1
   IF (k > maxitm) GO TO 30
   DO  j = k,maxitm
     hdr(j,i) = BLANK
   END DO
   30 nshft = nshft + nbpc
 END DO
 
 line  = nlpp + 1
 m0009 = 1023
 m1019 = lshift(1023,10)
 m2029 = lshift(1023,20)
 imask = lshift(1,30)
 
!     LOOP THROUGH DIT
 
 DO  jmkn = 1,ditsiz,2
   i = (jmkn-1)/2 + 1
   CALL fdit (i,k)
   ssname(1) = buf(k  )
   ssname(2) = buf(k+1)
   IF (ssname(1) == BLANK .AND. ssname(2) == BLANK) CYCLE
   CALL fmdi (i,k)
   
!     TEST TYPE BITS IN MDI
   
   DO  it = 2,ntype
     ibit = andf(buf(k+1),lshift(1,31-it))
     IF (ibit /= 0) GO TO 50
   END DO
   it = 1
   50 is = andf(buf(k+1),imask)
   im = BLANK
   IF (is /= 0) im = image
   ss = rshift(andf(buf(k+1),m1019),10)
   ps = andf(buf(k+1),m0009)
   ll = rshift(andf(buf(k+2),m2029),20)
   cs = rshift(andf(buf(k+2),m1019),10)
   hl = andf(buf(k+2),m0009)
   
!     LOOP THROUGH MDI ENTRY FOR THIS SUBSTRUCTURE DETERMINING THE
!     SIZE OF EACH EXISTING ITEM.
   
   DO  j = 1,nitm
     jj = j + ifrst - 1
     IF (buf(k+jj) == 0) GO TO 60
     inum = rshift(buf(k+jj),ihalf)*blksiz
     inum = ALOG10(FLOAT(inum)) + .3
     itm(j) = num(inum)
     IF (is /= 0 .AND. item(4,j) == 0) itm(j) = num(10)
     IF (ps /= 0 .AND. is == 0 .AND. item(5,j) == 0) itm(j) = num(10)
     CYCLE
     60 itm(j) = BLANK
   END DO
   
   line = line + 1
   IF (line <= nlpp) GO TO 90
   CALL page1
   line = line + 9 - 4
   WRITE  (nout,80) hdr
   80 FORMAT (//,26X,90HS u b s t r u c t u r e   o p e r a t i n g   f &
       &i l e   t a b l e   o f   c o n t e n t s , //,  &
       1H ,51X,26(a1,2X),a1,/1H ,51X,26(a1,2X),a1,/1H ,51X,26(a1,2X),a1,  &
       /,1H ,4X,12HSUBSTRUCTURE,35X,26(a1,2X),a1, /1H ,4X,3HNO.,3X,4HNAME, &
       4X,4HTYPE,3X,2HSS,3X,2HPS,3X,2HLL,3X,2HCS,3X,2HHL,4X,80(1H-)/)
   
   90 WRITE  (nout,100) i,ssname,im,TYPE(it),ss,ps,ll,cs,hl, (itm(l),l=1,nitm)
   100 FORMAT (2X,i6,2X,2A4,2X,a1,a2,5(1X,i4),4X,26(a1,2X),a1)
 END DO
 
!     PRINT SOF SPACE UTILIZATION MESSAGE
 
 line = line + 8
 IF (line > nlpp) CALL page1
 k    = sofsiz(k)
 nblk = 0
 DO  i = 1,nfiles
   nblk = nblk + filsiz(i)
 END DO
 iper = (avblks*100)/nblk
 WRITE  (nout,120) k,avblks,iper,hiblk
 120 FORMAT (//,51X,80HSIZE of item is given in powers of ten   (0 indi&
     &cates DATA is stored in primary) ,/,  &
     26H0*** unused SPACE on sof = ,i9,7H words.  ,/,  &
     22X,                   4HOR = ,i9,8H blocks. ,/,  &
     22X,                   4HOR = ,i9,9H percent.,/,  &
     26H0*** highest BLOCK used  = ,i9)
 line = nlpp
 RETURN
END SUBROUTINE softoc
