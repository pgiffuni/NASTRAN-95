SUBROUTINE eqsout
     
!     THIS ROUTINE WRITES THE CONNECTION TRACE FOR A NEWLY COMBINED
!     PSEUDOSTRUCTURE.
 
 EXTERNAL        lshift,rshift,andf,orf
 INTEGER :: z,score,orf,nbot(7),ntop(7),cnam
 INTEGER :: combo,iords(2)
 INTEGER :: words(6),ihd(64),string(32),andf,rshift
 COMMON /cmb002/ junk(5),score
 COMMON /cmb003/ combo(7,5),conset,iauto,toler,npsub,conect,tran,  &
     mcon,restct(7,7),isort,origin(7,3),iprint
 COMMON /cmb004/ tdat(6),nipnew,cnam(2)
 COMMON /zzzzzz/ z(1)
 COMMON /output/ ititl(96),ihead(96)
 COMMON /machin/ mach,ihalf
 DATA ihd/ 10*4H     , 4H   s , 4HUMMA , 4HRY o , 4HF ps , 4HEUDO ,  &
       4HSTRU , 4HCTUR , 4HE co , 4HNNEC , 4HTIVI ,  &
       4HTIES , 12*4H           , 4HINTE , 4HRNAL , 4H   i , 4HNTER ,  &
       4HNAL  , 4H  de , 4HGREE , 4HS of , 4H  ** ,5*4H****, 4H p s , 4H e u ,  &
       4H d o , 4H s t , 4H r u , 4H c t , 4H u r ,  &
       4H e   , 4H n a , 4H m e , 4H s * ,3*4H****, 3*4H      /
   DATA words / 4HPOIN , 4HT no , 4HFREE , 4HDOM  , 4HDOF  , 4HNO   /
   DATA iblank,nheqss  / 4H     , 4HEQSS /
   
   IF (andf(rshift(iprint,11),1) /= 1) RETURN
   CALL sfetch (cnam,nheqss,1,itest)
   CALL suread (z(score),-1,nout,itest)
   ncomp = z(score+2)
   nwrd  = nout - 4
   isteqs= score + nwrd
   
!     MOVE COMPONENT SUBSTRUCTURE NAMES INTO FIRST NWRD OF OPEN CORE.
   
   DO  i=1,nwrd
     ii = i - 1
     z(score+ii) = z(score+ii+4)
   END DO
   DO  i=1,32
     string(i) = iblank
   END DO
   CALL push (words(1),string, 5,8,0)
   CALL push (words(5),string,17,8,0)
   CALL push (words(3),string,29,8,0)
   DO  i=1,npsub
     iords(1) = combo(i,1)
     iords(2) = combo(i,2)
     loc = 39+11*(i-1)
     CALL push (iords(1),string,loc,8,0)
   END DO
   DO  i=1,64
     ihead(i) = ihd(i)
   END DO
   DO  i=65,96
     ihead(i) = string(i-64)
   END DO
   CALL page
   
!     COMPUTE FIRST AND LAST COMPONENT SUBSTRUCTURE ID NUMBERS
!     FOR EACH PSEUDOSTRUCTURE.
   
   nbot(1) = 1
   DO  i=1,npsub
     ntop(i) = nbot(i) + combo(i,5) - 1
     ii = i + 1
     IF (i == npsub) CYCLE
     nbot(ii) = ntop(i) + 1
   END DO
   
!     READ EQSS INTO OPEN CORE STARTING AT LOCATION ISTEQS
   
   jj    = 0
   icomp = 0
   180 icomp = icomp + 1
   IF (icomp > ncomp) GO TO 140
   170 CALL suread (z(isteqs+jj+1),3,nout,itest)
   SELECT CASE ( itest )
     CASE (    1)
       GO TO 130
     CASE (    2)
       GO TO 120
     CASE (    3)
       GO TO 140
   END SELECT
   130 CONTINUE
   
!     NORMAL ROUTE - PROCESS ENTRIES
   
   z(isteqs+jj) = icomp
   DO  j=1,npsub
     IF (icomp >= nbot(j) .AND. icomp <= ntop(j)) EXIT
   END DO
   150 z(isteqs+jj) = orf(lshift(j,ihalf),z(isteqs+jj))
   jj = jj + 4
   GO TO 170
   120 GO TO 180
   
!     SORT ON INTERNAL POINT NUMBER
   
   140 CONTINUE
   z(isteqs+jj  ) = 0
   z(isteqs+jj+1) = 0
   z(isteqs+jj+2) = 0
   z(isteqs+jj+3) = 0
   CALL sort (0,0,4,3,z(isteqs),jj)
   ii   = 1
   isil = 1
   DO  i=1,jj,4
     IF (z(isteqs+i+1) /= z(isteqs+i+5)) GO TO 210
     ii = ii + 1
     CYCLE
     210 iw = 4*ii
     ioff = i - 1 - 4*(ii-1)
     icode = z(isteqs+ioff+3)
     CALL decode (icode,string,ndof)
     CALL eqout1 (z(isteqs+ioff),iw,z(score),nwrd,isil)
     isil = isil + ndof
     ii   = 1
   END DO
   RETURN
 END SUBROUTINE eqsout
