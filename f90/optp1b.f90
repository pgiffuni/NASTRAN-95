SUBROUTINE optp1b (elt,elop,ele,pr)
     
 
 INTEGER, INTENT(IN OUT)                  :: elt(1)
 INTEGER, INTENT(IN OUT)                  :: elop(2,1)
 INTEGER, INTENT(IN OUT)                  :: ele(1)
 INTEGER, INTENT(OUT)                     :: pr(1)
 INTEGER :: count,ect,sysbuf,  &
     outtap,ycor,prcor,prc,NAME(2),card(2),elcr,elpt, pid,prpt,prpt1,b1p1
 COMMON /BLANK / skp1(2),count,skp2(2),ycor,b1p1,npow,  &
     nelw,nwdse,nprw,nwdsp,skp3, skp4(2),ect,skp5(4),numelm,itype(21)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /zzzzzz/ x(1)
 COMMON /optpw1/ prcor,prc(2)
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 COMMON /system/ sysbuf,outtap
 COMMON /names / nrd,noeor,nwrt,nweor
 DATA    NAME  / 4H opt,4HP1B  /
 
 
 ieop  = 1
 ides  = elop(1,ieop  )
 idee  = elop(1,ieop+1)
 prpt  = 1
 prpt1 = 1
 elop(2,1) = 1
 
!     IN CASE OF ERROR SET PRC(1)
 
 prc(1) = -1
 
 DO  k = 1,numelm
   nele = (idee-ides)/nwdse
   IF (nele < 0) THEN
     GO TO    30
   ELSE IF (nele == 0) THEN
     GO TO   120
   END IF
   
   10 idx = incr*(itype(k)-1)
   idp = idx + 4
   card(1) = NE(idp  )
   card(2) = NE(idp+1)
   IF (NE(idp+2) > prcor) GO TO 150
   CALL locate (*160,x(b1p1),card(1),i)
   
!     SEQUENTIAL ELEMENT SEARCH
   
   npr  = 0
   elpt = ides
   elcr = ele(elpt)
   
   20 CALL READ (*160,*160,ect,prc,NE(idp+2),noeor,i)
   IF (prc(1)-elcr < 0.0) THEN
     GO TO    20
   ELSE IF (prc(1)-elcr == 0.0) THEN
     GO TO    50
   END IF
   
!     LOGIC OR FILE FAILURE
   
   30 CALL page2 (-2)
   WRITE  (outtap,40) sfm,itype(k),prc(1),NAME
   40 FORMAT (a25,' 2297, INCORRECT LOGIC FOR ELEMENT TYPE',i4,  &
       ', ELEMENT',i8,2H (,2A4,2H).)
   GO TO 170
   
!     ELEMENT ID IN CORE .EQ. ECT ID - ELEMENT TO BE OPTIMIZED
   
   50 pid = prc(2)
   card(1) = pid
   card(2) = ele(elpt+4)
   
!     TEST FOR CORE NEEDED AFTER EXPANDING TO NWDSP WORDS
   
   IF (prpt1+nwdsp*(npr/2+1) > ycor) GO TO 180
   CALL bishel (*60,card,npr,2,pr(prpt1))
   60 ele(elpt+4) = pid
   elpt = elpt + nwdse
   IF (elpt >= idee) GO TO 70
   elcr = ele(elpt)
   GO TO 20
   
!     NEW ELEMENT TYPE COMING
   
   70 CALL fread (ect,0,0,nweor)
   
!     EXPAND PROPERTIES TO NWDSP WORDS/PROPERTY
   
   nx = npr/2
   IF (nx-1 < 0) THEN
     GO TO    30
   ELSE IF (nx-1 == 0) THEN
     GO TO   100
   END IF
   80 CONTINUE
   DO  i = 1,nx
     j = nx - i
     l = prpt1 + j*nwdsp
     m = prpt1 + j*2
     pr(l  ) = pr(m  )
     pr(l+1) = pr(m+1)
   END DO
   
   100 prpt = prpt1 + nx*nwdsp
   
!     PLACE POINTERS IN ELEMENT ARRAY
   
   l = idee - 1
   DO  i = ides,l,nwdse
     kid = ele(i+4)
     CALL bisloc (*30,kid,pr(prpt1),nwdsp,nx,j)
     ele(i+4) = j
   END DO
   
!     SETUP FOR NEXT ELEMENT
   
   120 ieop = ieop + 1
   elop(2,ieop) = prpt
   prpt1 = prpt
   ides  = idee
   IF (ieop > npow) GO TO 140
   idee = elop(1,ieop+1)
 END DO
 
 
 140 nprw = prpt - 1
 RETURN
 
!     ERRORS
 
!     INSUFFICIENT CORE IN /OPTPW1/ OR /XXOPT1/
 
 150 count = -1
 GO TO 140
 
!     FILE ERRORS
 
 160 CALL mesage (-7,ect,NAME)
 170 prpt = 1
 GO TO 140
 
!     INSUFFICIENT CORE
 
 180 CALL page2 (-2)
 WRITE  (outtap,190) ufm,NAME,b1p1,pid
 190 FORMAT (a23,' 2298, INSUFFICIENT CORE ',2A4,1H(,i10,'), PROPERTY', i9)
 GO TO 150
END SUBROUTINE optp1b
