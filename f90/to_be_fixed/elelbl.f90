SUBROUTINE elelbl (gplst,x,u,deform,buf1)
     
 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(IN OUT)                  :: buf1
 LOGICAL :: solid
 INTEGER :: BLANK    , br       ,q4       ,t3      ,  &
     chr      ,crew     , hb       ,twod    ,  &
     ect      ,elid     ,elidp(2) ,elsets   ,eltype  ,  &
     gpts     , lbl(10)  ,lblp(8)  ,rdrew   ,  &
     pid      ,plabel   ,pltflg   ,pset     ,offset
 REAL :: infnty   ,LEN      ,ma       ,maxlen   ,mb      , minslp
 COMMON /BLANK / skp(3),pltflg,skp1(6),skp2(2),elsets,casecc(5),ect
 COMMON /system/ ksystm(40),ncpw
 COMMON /names / norew,rdrew,skpn(2),crew
 COMMON /gpta1 / ntyps,last,incr,NE(1)
 COMMON /char94/ chr(60)
 COMMON /pltdat/ skpplt(20),skpa(3),cntchr(2)
 COMMON /drwdat/ pset,plabel
 COMMON /pltscr/ ncor,xy(2,8),gpts(4)
 
 DATA   BLANK  , infnty,slpmax   / 1H ,1.e3,5. /, pid / 4      /,  &
     itetra / 2HTE   /, iect  / 4HECT   / , hb     / 2HHB   /,  &
     br     / 2HBR   /, q4    / 2HQ4    /,  t3     / 2HT3   /
 
 np = 0
 CALL tipe (0,0,0,0,0,-1)
 cntx = cntchr(1)
 cnty = cntchr(1) + (cntchr(2) - cntchr(1))/2.
 
!   . CHECK IF PROPERTY ID IS TO BE TYPED NEXT TO ELEMENT ID
 
 IF (plabel /= pid) GO TO 40
 IF (pltflg <   0) GO TO 10
 CALL preloc (*10,gplst(buf1),ect)
 CALL fname  (ect,gpts(1))
 IF (gpts(1) == iect) GO TO 20
 CALL CLOSE (ect,crew)
 10 plabel = pid - 1
 GO TO 40
 20 CALL delset
 lpid = 0
 
!   . READ THE ELEMENT TYPE + NUMBER OF GRID POINTS / ELEMENT OF THIS
!     TYPE.
 
 30 IF (lpid > 0) CALL fwdrec(*40,ect)
 40 CALL READ (*200,*200,elsets,eltype,1,0,i)
 CALL fread (elsets,ngpel,1,0)
 twod = 0
 IF (ngpel > 2) twod = 1
 ngpel = IABS (ngpel)
 solid =.false.
 IF ((eltype == itetra .OR. ngpel > 4) .AND. eltype /= hb) solid = .true.
!-----
!   . REJECT ELEMENTS WITH 0 OR MORE THAN --NCOR-16-- GRID POINTS
 
 IF (ngpel > 1 .AND. ngpel < ncor-13) GO TO 60
 50 CALL fread (elsets,elid,1,0)
 IF (elid <= 0) GO TO 40
 CALL fread (elsets,0,-1,0)
 CALL fread (elsets,0,-ngpel,0)
 GO TO 50
 60 CONTINUE
!-----
 IF (plabel /= pid) GO TO 90
 j = 16
 DO  i = 1,ntyps
   IF (NE(j) == eltype) GO TO 80
   j = j + incr
 END DO
 GO TO 90
 80 lpid = j - 12
 IF (NE(lpid+2) <= 0) GO TO 90
 npid = NE(lpid+2)
 CALL locate (*90,gplst(buf1),NE(lpid),gpts(1))
 GO TO 100
 90 lpid = 0
 
 100 ngpel1 = ngpel + 1
 offset = 0
 IF (eltype == br) offset = 6
 IF (eltype == q4 .OR. eltype == t3) offset = 1
 
!     READ AN ELEMENT ID + ITS GRID POINTS.
 
 102 CALL fread (elsets,elid,1,0)
 IF (eltype == hb) ngpel = 8
 IF (elid   <=  0) GO TO 30
 CALL fread (elsets,0,-1,0)
 CALL fread (elsets,gpts(1),ngpel,0)
 IF (offset > 0) CALL fread (elsets,0,-offset,0)
 IF (eltype /= hb) GO TO 1028
 DO  i = 2,4
   IF (gpts(i) == 0.0) THEN
     GO TO  1025
   ELSE
     GO TO  1023
   END IF
 END DO
 i = 5
 1025 ngpel = i - 1
 1028 CONTINUE
 k = elid
 nl = 0
 DO  i = 1,8
   j = elid/10**(8-i)
   IF (j == 0 .AND. nl == 0) CYCLE
   nl = nl + 1
   lbl(nl) = chr(j+1)
   elid = elid - j*10**(8-i)
 END DO
 lbl(nl+1) = khrfn1(BLANK,1,eltype,1)
 lbl(nl+2) = khrfn1(BLANK,1,eltype,2)
 nl = nl + 2
 
!   . DECODE PROPERTY ID
 
 IF (lpid <= 0) GO TO 105
 1040 CALL READ (*1041,*1041,ect,elidp,2,0,i)
 CALL fread (ect,0,-(npid-2),0)
 IF (elidp(1) == k) GO TO 1042
 GO TO 1040
 1041 lpid = -1
 GO TO 105
 
!   . ELEMENT PROPERTY FOUND
 
 1042 k  = 10000000
 np = 0
 DO  i = 1,8
   j = elidp(2)/k
   IF (j == 0 .AND. np == 0) GO TO 1043
   np = np + 1
   lblp(np) = chr(j+1)
   elidp(2) = elidp(2) - j*k
   1043 k = k/10
 END DO
 
 105 CONTINUE
 
!   . SET UP THE COORDINATES OF THE GRID POINTS
 
 DO  i = 1,ngpel
   j = gpts(i)
   j = IABS(gplst(j))
   IF (deform /= 0) GO TO 106
   xx = x(2,j)
   yy = x(3,j)
   GO TO 107
   106 xx = u(1,j)
   yy = u(2,j)
   107 IF (solid) GO TO 1071
   xy(1,i) = xx
   xy(2,i) = yy
   j = ngpel + i
   xy(1,j) = xx
   xy(2,j) = yy
   CYCLE
   1071 IF (i > 2) GO TO 1072
   xy(1,i) = xx
   xy(2,i) = yy
   IF (i /= 1) GO TO 1072
   xy(1,3) = 0.0
   xy(2,3) = 0.0
   1072 xy(1,3) = xx + xy(1,3)
   xy(2,3) = yy + xy(2,3)
 END DO
 
 IF (solid) GO TO 160
 IF (twod  /= 0) GO TO 110
 IF (ngpel == 2) GO TO 125
 k = 3
 GO TO 120
 
!     FIND THE BASE OF THIS POLYGON = LONGEST SIDE (IF MORE THAN ONE
!     LONGEST SIDE, CHOOSE FROM THEM THE SIDE OF SMALLEST SLOPE).
 
 110 maxlen = 0.
 DO  i = 1,ngpel
   xx = xy(1,i+1) - xy(1,i)
   yy = xy(2,i+1) - xy(2,i)
   LEN = xx**2 + yy**2
   IF (xx /= 0.) GO TO 111
   slp = infnty
   GO TO 112
   111 slp = ABS(yy/xx)
   112 IF (maxlen-LEN < 0) THEN
     GO TO   113
   ELSE IF (maxlen-LEN == 0) THEN
     GO TO   114
   ELSE
     GO TO   116
   END IF
   113 maxlen = LEN
   GO TO 115
   114 IF (slp >= minslp) CYCLE
   115 k = i
   minslp = slp
 END DO
 
 IF (k == 1) GO TO 122
 120 DO  i = 1,ngpel1
   xy(1,i) = xy(1,k)
   xy(2,i) = xy(2,k)
   k = k + 1
 END DO
 122 IF (ngpel == 6) GO TO 140
 IF (ngpel-3 < 0) THEN
   GO TO   125
 ELSE IF (ngpel-3 == 0) THEN
   GO TO   140
 ELSE
   GO TO   150
 END IF
 
!     LINE ELEMENT.
 
 125 xx = xy(1,2) - xy(1,1)
 IF (xx == 0.) GO TO 126
 yy = xy(2,2) - xy(2,1)
 slp= yy/xx
 GO TO 127
 126 slp= infnty
 127 xc = (xy(1,1) + xy(1,2))/2.
 yc = (xy(2,1) + xy(2,2))/2.
 
 IF (ABS(slp)-1. > 0.0) THEN
   GO TO   129
 END IF
 128 yc = yc + cnty
 GO TO 175
 129 IF (ABS(slp)-slpmax < 0.0) THEN
   GO TO   130
 ELSE
   GO TO   131
 END IF
 130 xc = xc - SIGN(cntx,slp)
 GO TO 175
 131 xc = xc + cntx
 GO TO 175
 
!     TRIANGULAR ELEMENT.  POINTS 1+2 ARE THE BASE - POINT 3 THE APEX.
 
 140 xc = (xy(1,1) + xy(1,2) + xy(1,3))/3.
 yc = (xy(2,1) + xy(2,2) + xy(2,3))/3.
 GO TO 175
 
!     QUADRILATERAL ELEMENT.
 
 150 xx = (xy(1,3)+xy(1,4)) - (xy(1,1)+xy(1,2))
 IF (xx /= 0.) GO TO 151
 ma = infnty
 GO TO 152
 151 yy = (xy(2,3)+xy(2,4)) - (xy(2,1)+xy(2,2))
 ma = yy/xx
 ba = (xy(2,1)+xy(2,2))/2. - ma*(xy(1,1)+xy(1,2))/2.
 152 xx = (xy(1,2)+xy(1,3)) - (xy(1,1)+xy(1,4))
 IF (xx /= 0.) GO TO 153
 mb = infnty
 GO TO 155
 153 yy = (xy(2,2)+xy(2,3)) - (xy(2,1)+xy(2,4))
 mb = yy/xx
 bb = (xy(2,1)+xy(2,4))/2. - mb*(xy(1,1)+xy(1,4))/2.
 
 155 IF (ABS(ma) >= infnty) GO TO 156
 IF (ABS(mb) >= infnty) GO TO 157
 IF (mb == ma) GO TO 158
 xc = (ba-bb)/(mb-ma)
 yc = ma*xc + ba
 GO TO 175
 156 xc = (xy(1,1) + xy(1,2))/2.
 yc = mb*xc + bb
 GO TO 175
 157 xc = (xy(1,1) + xy(1,4))/2.
 yc = ma*xc + ba
 GO TO 175
 158 xc = (xy(1,3) + xy(1,4) + xy(1,2)+xy(1,1))/4.0
 yc = (xy(2,3) + xy(2,4) + xy(2,2)+xy(2,1))/4.0
 GO TO 175
 
!   . ELEMENTS WITH MORE THAN FOUR GRIDS
 
 160 xc = xy(1,3)/FLOAT(ngpel)
 yc = xy(2,3)/FLOAT(ngpel)
 GO TO 175
 
!     SETUP THE STRAIGHT LINE EQUATION OF THE LINE ON WHICH THE ELEMENT
!     LABEL IS TO BE TYPED - Y=MX+B.
 
 175 xx = xy(1,2) - xy(1,1)
 IF (xx == 0.) GO TO 176
 yy = xy(2,2) - xy(2,1)
 slp= yy/xx
 b  = yc - xc*slp
 GO TO 180
 176 slp= infnty
 
!     TYPE THE ELEMENT LABEL (NL CHARACTERS)
 
 180 zz = nl/2
 IF (nl/2 == (nl+1)/2) zz = zz - .5
 absslp = ABS(slp)
 cc = cntx
 IF (absslp >= slpmax) cc = cnty
 k = MAX0(nl,np)
 
 DO  i = 1,k
   xx = cc*(zz - FLOAT(i-1))
   IF (absslp > 1.) GO TO 181
   xx = xc - xx
   yy = slp*xx + b
   GO TO 190
   181 IF (absslp >= slpmax) GO TO 182
   yy = SIGN(1.,slp)
   GO TO 183
   182 yy = -1.
   183 yy = yc - yy*xx
   IF (absslp >= infnty) GO TO 184
   xx = (yy-b)/slp
   GO TO 190
   184 xx = xc
   
   
!     OFFSET THE HB LABEL AND PROPERTY ID IF ANY WHEN TIPE LABEL
   
   190 IF (eltype /= hb) GO TO 1905
   jtj = 2
   IF (absslp < slpmax) yy = yy - jtj*cc
   IF (absslp >= slpmax) xx = xx + jtj*cc
   1905 IF (nl >= i) CALL tipe (xx,yy,1,lbl(i),1,0)
   IF (lpid <= 0) CYCLE
   IF (np   < i) CYCLE
   IF (absslp < slpmax) yy = yy - 2.*cc
   IF (absslp >= slpmax) xx = xx + 2.*cc
   CALL tipe (xx,yy,1,lblp(i),1,0)
 END DO
 GO TO 102
 
 200 CALL tipe (0,0,0,0,0,1)
 IF (plabel == pid) CALL CLOSE (ect,crew)
 RETURN
END SUBROUTINE elelbl
