SUBROUTINE pthbdy
     
!     PTHBDY MODIFIES THE SIL,ECT,EQEXIN AND BGBDT FOR CHBDY ELEMENTS
!     SO THEY CAN BE PLOTTED.
 
!     THE SIL BGPDT AND EQEXIN OUTPUT LOOKALIKES ARE ADD ON FILES
!     THE ECT OUTPUT FILE HAS THE CHBDY FLAG SET NEGATIVE
!     SO PLTSET CAN TELL THE ECTS APART IT ALSO HAS THE NEW GRID POINTS
 
 INTEGER :: NAME(2),ept,ect,sil,ieq,bgpdt,hect,hsil,oeq,  &
     hbgpdt,geom2,scr1,scr2,flag,iz(1),sysbuf,out,  &
     cbs(20),osil,file1,file2,view(2),chbdy(2),  &
     phbdy(2),buf1,buf2,buf3,buf4,buf5,buf6,trl(7)
 DIMENSION       nsil(7),neq(14),tem(3),e(3),v(3),r21(3)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / nhbdy,mesh(2)
 COMMON /system/ sysbuf,out,dum(6),nlpp
 COMMON /condas/ pi
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (iz(1),z(1))
 DATA    geom2 , ect, ept, sil, ieq, bgpdt                   /  &
     101   , 102, 103, 104, 105,   106                   /
 DATA    hect  , hsil, oeq, hbgpdt, scr1, scr2               /  &
     201   , 202,  203,    204,  301,  302               /
 DATA    iyes  , no    , NAME  ,         nphbdy, nview, ncb2 /  &
     4HYES , 4HNO  , 4HPLNB, 4HDY  , 7     , 6    , 15   /
 DATA    view  ,         chbdy  ,        phbdy        , nect /  &
     2606  , 26    , 4208,42,        2502,25      , 15   /
 
!     PRINT FLAG CHBDY FLAG
 
 iprt  = 0
 IF (mesh(1) == iyes) iprt = 1
 nhbdy = -1
 line  = nlpp
 
!     INITIALIZE
 
 buf1 = korsz(z(1)) - sysbuf
 buf2 = buf1 - sysbuf  - 1
 buf3 = buf2 - sysbuf
 buf4 = buf3 - sysbuf
 buf5 = buf4 - sysbuf
 buf6 = buf5 - sysbuf
 CALL preloc (*1000,z(buf1),geom2)
 CALL locate (*1000,z(buf1),chbdy,n)
 
!     MAKE A SCRATCH FILE WITH EID AF DISLIN FOR CHBDY
 
 ipn  = 0
 ivew = 0
 CALL preloc (*25,z(buf2),ept)
 file1 = ept
 CALL locate (*15,z(buf2),phbdy,n)
 CALL READ (*1002,*10,ept,z(1),buf3,0,n)
 GO TO 1008
 10 ipn = n
 15 ipv = ipn + 1
 nrd = buf3 - ipv
 CALL locate (*25,z(buf2),view,n)
 CALL READ (*1002,*20,ept,z(ipv),nrd,0,n)
 GO TO 1008
 20 ivew = n
 25 CALL CLOSE (ept,1)
 CALL gopen (scr1,z(buf2),1)
 file1 = geom2
 30 CALL READ (*1002,*70,geom2,cbs,ncb2,0,n)
 tem(1) = 0.0
 tem(2) = 0.0
 IF (ipn == 0) GO TO 40
 DO  i = 1,ipn,nphbdy
   IF (cbs(2) /= iz(i)) CYCLE
   tem(1) = z(i+2)
   EXIT
 END DO
 40 IF (ivew    == 0) GO TO 65
 IF (cbs(15) == 0) GO TO 65
 DO  i = 1,ivew,nview
   IF (cbs(15) /= iz(ipn+i)) CYCLE
   tem(2) = z(ipn+i+5)
   IF (iprt ==    0) EXIT
   IF (line < nlpp) GO TO 50
   line = 1
   CALL page1
   WRITE  (out,45)
   45 FORMAT (1H0,17X,5HIDENT,8X,4HBETA,7X,5HGAMMA,9X,3HCAN,6X,6HCAN be,  &
       /6X,5HCHBDY,6X,6HNUMBER,8X,4HMESH,8X,4HMESH,7X,5HSHADE,6X,  &
       6HSHADED,5X,7HDISLIN ,/)
   50 nb = iyes
   ns = iyes
   IF (iz(ipn+i+1) == 0) nb = no
   IF (iz(ipn+i+2) == 0) ns = no
   line = line+1
   WRITE (out,55) cbs(1),iz(ipn+i),iz(ipn+i+3),iz(ipn+i+4),nb,ns, tem(2)
   55 FORMAT (1H ,4(i10,2X),6X,a4,8X,a4,2X,1P,e10.4)
   EXIT
 END DO
 65 CALL WRITE (scr1,tem,2,0)
 GO TO 30
 70 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1 ,1)
 CALL CLOSE (geom2,1)
 CALL gopen (scr1,z(buf1),0)
 trl(1) = sil
 CALL rdtrl (trl)
 osil   = trl(3)
 trl(1) = bgpdt
 CALL rdtrl (trl)
 nin = trl(2)
 nrd = 4*trl(2)
 IF (5*sysbuf+nrd+50 > buf1) GO TO 1008
 
!     FIND CHBDY CARDS COPY ECT TO CHBDY CARDS
 
 CALL gopen (ect,z(buf2),0)
 CALL gopen (hect,z(buf3),1)
 file1 = ect
 80 CALL READ (*1000,*1000,ect,cbs,3,0,n)
 CALL WRITE (hect,cbs,3,0)
 IF (cbs(1) == chbdy(1) .AND. cbs(2) == chbdy(2)) GO TO 95
 
!     DUPE REST OF RECORD
 
 85 CALL READ (*1002,*90,ect,z(1),buf6-1,0,n)
 CALL WRITE (hect,z(1),buf6-1,0)
 GO TO 85
 90 CALL WRITE (hect,z(1),n,1)
 GO TO 80
 
!     COPY SIL EQEXIN TO NEW FILES
 
 95 icore = buf6 - 1
 ilft  = 1
 file1 = sil
 file2 = hsil
 n     = buf4
 leq   = 0
 100 CALL gopen (file1,z(buf6),0)
 CALL gopen (file2,z(n),1)
 105 CALL READ (*1002,*120,file1,z(ilft),icore,0,m)
 IF (file1 /= ieq) GO TO 115
 DO  i = 1,icore,2
   IF (iz(i) > leq) leq = iz(i)
 END DO
 115 CALL WRITE (file2,z(ilft),icore,0)
 GO TO 105
 120 IF (file1 /= ieq) GO TO 130
 DO  i = 1,m,2
   IF (iz(i) > leq) leq = iz(i)
 END DO
 130 CALL WRITE (file2,z(ilft),m,0)
 CALL CLOSE (file1,1)
 IF (n == buf5) GO TO 150
 file1 = ieq
 file2 = oeq
 n     = buf5
 GO TO 100
 
!     BRING IN BGPDT
 
 150 CALL gopen (bgpdt,z(buf6),0)
 file1 = bgpdt
 CALL READ (*1002,*1002,bgpdt,z(1),nrd,0,n)
 CALL CLOSE (bgpdt,1)
 
!     FINALLY TIME TO GO TO WORK
 
 nhbdy = 0
 nngp  = 0
 nbgp  = nrd + 1
 DO  i = 1,24
   iz(nrd+i) = 0
 END DO
 file1 = ect
 CALL gopen (scr2,z(buf6),1)
 160 CALL READ (*1002,*350,ect,cbs,nect,0,n)
 cbs(nect) = 0.0
 IF (cbs(3) > 6) cbs(3) = 3
 CALL READ (*1002,*1002,scr1,tem,2,0,n)
 flag  = cbs(3)
 nhbdy = nhbdy + 1
 SELECT CASE ( flag )
   CASE (    1)
     GO TO 200
   CASE (    2)
     GO TO 250
   CASE (    3)
     GO TO 260
   CASE (    4)
     GO TO 270
   CASE (    5)
     GO TO 280
   CASE (    6)
     GO TO 260
 END SELECT
 
!     POINT
 
 
!     BGPDT DATA FOR POINT
 
 200 i1 = (cbs(4)-1)*4 + 2
 itry = 1
 e(1) = 0.0
 e(2) = 0.0
 e(3) = 0.0
 CALL sapb (cbs(12),e,v)
 CALL sanorm (*400,v)
 e(1) = 1.0
 205 xl = sadotb(v,e)
 r21(1) = e(1) - xl*v(1)
 r21(2) = e(2) - xl*v(2)
 r21(3) = e(3) - xl*v(3)
 xl = sadotb(r21,r21)
 IF (xl  > .2) GO TO 210
 IF (itry == 2) GO TO 400
 itry = 2
 e(1) = 0.0
 e(2) = 1.0
 GO TO 205
 210 CALL sanorm (*215,r21)
 215 CALL saxb (v,r21,e)
 xl = 0.0
 IF (tem(1) /= 0.0) xl = SQRT(tem(1)/pi)
 zs3 = .8660254
 z(nbgp+ 1) = z(i1  ) + xl*r21(1)
 z(nbgp+ 2) = z(i1+1) + xl*r21(2)
 z(nbgp+ 3) = z(i1+2) + xl*r21(3)
 z(nbgp+ 5) = z(i1  ) + xl*( .5*r21(1) + zs3*e(1))
 z(nbgp+ 6) = z(i1+1) + xl*( .5*r21(2) + zs3*e(2))
 z(nbgp+ 7) = z(i1+2) + xl*( .5*r21(3) + zs3*e(3))
 z(nbgp+ 9) = z(i1  ) + xl*(-.5*r21(1) + zs3*e(1))
 z(nbgp+10) = z(i1+1) + xl*(-.5*r21(2) + zs3*e(2))
 z(nbgp+11) = z(i1+2) + xl*(-.5*r21(3) + zs3*e(3))
 z(nbgp+14) = z(i1+1) - xl* r21(2)
 z(nbgp+15) = z(i1+2) - xl* r21(3)
 z(nbgp+17) = z(i1  ) + xl*(-.5*r21(1) - zs3*e(1))
 z(nbgp+18) = z(i1+1) + xl*(-.5*r21(2) - zs3*e(2))
 z(nbgp+19) = z(i1+2) + xl*(-.5*r21(3) - zs3*e(3))
 z(nbgp+21) = z(i1  ) + xl*(+.5*r21(1) - zs3*e(1))
 z(nbgp+22) = z(i1+1) + xl*(+.5*r21(2) - zs3*e(2))
 z(nbgp+23) = z(i1+2) + xl*(+.5*r21(3) - zs3*e(3))
 z(nbgp+25) = z(i1  ) + xl*v(1)
 z(nbgp+26) = z(i1+1) + xl*v(2)
 z(nbgp+27) = z(i1+2) + xl*v(3)
 nngp= nngp + 7
 ns  = 7
 nea = 14
 nb  = 28
 m   = 7
 n   = 5
 220 nn  = 1
 DO  i = 1,m
   leq = leq + 1
   nin = nin + 1
   nsil(i) = nin
   neq(nn) = leq
   neq(nn+1) = nin
   cbs(n   ) = nin
   nn  = nn + 2
   n   = n  + 1
 END DO
 GO TO 300
 
!     LINE
 
 
!     BGPDT DATA FOR LINE
 
 250 i1 = (cbs(4)-1)*4 + 2
 i2 = (cbs(5)-1)*4 + 2
 CALL samb (z(i2),z(i1),r21)
 xl = sadotb(r21,r21)
 IF (xl == 0.0) GO TO 400
 x1 = sadotb(r21,cbs(12))
 xl = x1/xl
 e(1) = xl*r21(1)
 e(2) = xl*r21(2)
 e(3) = xl*r21(3)
 CALL samb (cbs(12),e,v)
 CALL sanorm (*400,v)
 CALL saxb   (v,r21,e)
 CALL sanorm (*400,e)
 d  = tem(2)
 af = tem(1)*.5
 z(nbgp+ 1) = z(i1  ) + d*v(1) - af*e(1)
 z(nbgp+ 2) = z(i1+1) + d*v(2) - af*e(2)
 z(nbgp+ 3) = z(i1+2) + d*v(3) - af*e(3)
 z(nbgp+ 5) = z(i2  ) + d*v(1) - af*e(1)
 z(nbgp+ 6) = z(i2+1) + d*v(2) - af*e(2)
 z(nbgp+ 7) = z(i2+2) + d*v(3) - af*e(3)
 z(nbgp+ 9) = z(i2  ) + d*v(1) + af*e(1)
 z(nbgp+10) = z(i2+1) + d*v(2) + af*e(2)
 z(nbgp+11) = z(i2+2) + d*v(3) + af*e(3)
 z(nbgp+13) = z(i1  ) + d*v(1) + af*e(1)
 z(nbgp+14) = z(i1+1) + d*v(2) + af*e(2)
 z(nbgp+15) = z(i1+2) + d*v(3) + af*e(3)
 z(nbgp+17) = z(i1  ) + d*v(1) + .5*r21(1)
 z(nbgp+18) = z(i1+1) + d*v(2) + .5*r21(2)
 z(nbgp+19) = z(i1+2) + d*v(3) + .5*r21(3)
 z(nbgp+21) = z(nbgp+17) + 2.*af*v(1)
 z(nbgp+22) = z(nbgp+18) + 2.*af*v(2)
 z(nbgp+23) = z(nbgp+19) + 2.*af*v(3)
 nngp= nngp + 6
 ns  = 6
 nea = 12
 nb  = 24
 m   = 6
 n   = 6
 GO TO 220
 
!     REV  OR  ELIP   DO NOTHING
 
 260 GO TO 310
 
!     AREA3
 
!     BGPDT DATA FOR AREA3
 
 270 i1 = (cbs(4)-1)*4 + 2
 i2 = (cbs(5)-1)*4 + 2
 i3 = (cbs(6)-1)*4 + 2
 CALL samb (z(i2),z(i1),e)
 CALL samb (z(i3),z(i1),v)
 CALL saxb (e,v,e)
 CALL sanorm (*400,e)
 CALL samb (z(i2),z(i1),v)
 x1 = sadotb(v,v)
 CALL samb (z(i3),z(i1),v)
 x2 = sadotb(v,v)
 CALL samb (z(i3),z(i2),v)
 x3 = sadotb(v,v)
 x1 = AMAX1(x1,x2)
 x1 = AMAX1(x1,x3)
 xl = .25* SQRT(x1)
 CALL sapb (z(i1),z(i2),v)
 CALL sapb (z(i3),v,v)
 z(nbgp+1) = v(1)/3.0
 z(nbgp+2) = v(2)/3.0
 z(nbgp+3) = v(3)/3.0
 z(nbgp+5) = z(nbgp+1) + xl*e(1)
 z(nbgp+6) = z(nbgp+2) + xl*e(2)
 z(nbgp+7) = z(nbgp+3) + xl*e(3)
 275 nngp= nngp + 2
 ns  = 2
 nea = 4
 nb  = 8
 leq = leq + 1
 nin = nin + 1
 n   = 7
 IF (flag == 5) n = 8
 nsil(1)= nin
 neq(1) = leq
 neq(2) = nin
 cbs(n) = nin
 leq = leq + 1
 nin = nin + 1
 nsil( 2) = nin
 neq ( 3) = leq
 neq ( 4) = nin
 cbs(n+1) = nin
 cbs(n+2) = nin
 cbs(n+3) = nin
 IF (flag == 4) cbs(n+4) = nin
 GO TO 300
 
!     AREA4
 
!     BGPDT DATA FOR AREA4
 
 280 i1 = (cbs(4)-1)*4 + 2
 i2 = (cbs(5)-1)*4 + 2
 i3 = (cbs(6)-1)*4 + 2
 i4 = (cbs(7)-1)*4 + 2
 CALL samb (z(i3),z(i1),e)
 CALL samb (z(i4),z(i2),v)
 CALL saxb (e,v,e)
 CALL sanorm (*400,e)
 CALL samb (z(i2),z(i1),v)
 x1 = sadotb(v,v)
 CALL samb (z(i3),z(i2),v)
 x2 = sadotb(v,v)
 CALL samb (z(i4),z(i3),v)
 x3 = sadotb(v,v)
 CALL samb (z(i4),z(i1),v)
 x4 = sadotb(v,v)
 x1 = AMAX1(x1,x2)
 x1 = AMAX1(x1,x3)
 x1 = AMAX1(x1,x4)
 xl = .25* SQRT(x1)
 CALL sapb (z(i1),z(i2),v)
 CALL sapb (v,z(i3),v)
 CALL sapb (v,z(i4),v)
 z(nbgp+1) = .25*v(1)
 z(nbgp+2) = .25*v(2)
 z(nbgp+3) = .25*v(3)
 z(nbgp+5) = z(nbgp+1) + xl*e(1)
 z(nbgp+6) = z(nbgp+2) + xl*e(2)
 z(nbgp+7) = z(nbgp+3) + xl*e(3)
 GO TO 275
 
!     ADD TO HSIL HEQEXIN  HECT
!     BGPDT
 
 300 CALL WRITE (hsil,nsil,ns,0)
 CALL WRITE (oeq,neq,nea,0)
 CALL WRITE (scr2,z(nbgp),nb,0)
 310 cbs(3) = -cbs(3)
 CALL WRITE (hect,cbs,nect,0)
 GO TO 160
 
!     END CLOSE FILES, WRITE NBGPDT, WRITE TRAILERS THEN FINISH ECT COPY
 
 350 CALL WRITE (hsil,0,0,1)
 CALL WRITE (oeq ,0,0,1)
 CALL WRITE (hect,0,0,1)
 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr2,1)
 CALL CLOSE (scr1,1)
 CALL CLOSE (hsil,1)
 CALL CLOSE (oeq ,1)
 CALL gopen (hbgpdt,z(buf1),1)
 CALL WRITE (hbgpdt,z,nrd,0)
 IF (nngp == 0) GO TO 380
 file1 = scr2
 CALL gopen (scr2,z(buf6),0)
 360 CALL READ (*1002,*370,scr2,z(1),buf6-1,0,n)
 CALL WRITE (hbgpdt,z(1),buf6-1,0)
 GO TO 360
 370 CALL WRITE (hbgpdt,z(1),n,1)
 CALL CLOSE (scr2,1)
 380 CALL CLOSE (hbgpdt,1)
 trl(1) = hbgpdt
 trl(2) = nrd/4  + nngp
 CALL wrttrl (trl)
 trl(1) = oeq
 CALL wrttrl (trl)
 trl(1) = hsil
 trl(3) = nngp + osil
 CALL wrttrl (trl)
 trl(1) = ect
 CALL rdtrl (trl)
 trl(1) = hect
 CALL wrttrl (trl)
 file1  = ect
 GO TO 80
 
!     BAD GEOMETRY FOR ELEMENT
 
 400 cbs(3) = -cbs(3)
 nhbdy  = nhbdy - 1
 WRITE  (out,410) uwm,cbs(1)
 410 FORMAT (a25,', CHBDY ELEMENT',i9,' HAS NO NORMAL OR BAD GEOMETRY',  &
     ' WHICH MAKES IT UNPLOTTABLE')
 GO TO 310
 
!     RETURN OR ERROR MESSAGES
 
 1000 CALL CLOSE (ect,1)
 CALL CLOSE (hect,1)
 CALL CLOSE (geom2,1)
 RETURN
 
 1002 CALL mesage (-2,0,file1)
 1008 CALL mesage (-8,0,NAME)
 RETURN
END SUBROUTINE pthbdy
