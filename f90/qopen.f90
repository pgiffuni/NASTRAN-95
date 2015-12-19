SUBROUTINE qopen ( *, namfil, buff, iop )
 
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
  
 INTEGER, INTENT(IN)                      :: namfil
 INTEGER, INTENT(IN)                      :: buff(10)
 INTEGER, INTENT(IN)                      :: iop
 COMMON / system / isysbf, dum1(77), idiag, dum2(21)
 INTEGER :: dname(2), itrl(7)
 DATA     init   / 0 /
 
 NAME   = namfil
 iocode = iop
 IF ( init /= 0 ) GO TO 10
 ibasbf = locfx( ibase )
 CALL dsiodd
 nbuff  = isysbf - 4
 nbfz   = 1
 IF ( lenwpb /= 0 ) nbfz   = nbuff / lenwpb + .1
 init   = 1
 10   IF ( IAND( idiag, 2**14 ) /= 0 ) CALL dsmsg ( 1 )
 locbuf = locfx( buff )
 indbas = locbuf - ibasbf + 1
 IF ( MOD( indbas,2 ) == 0 ) indbas = indbas + 1
 IF ( fcb( 2, ifilex ) == 0 ) GO TO 20
 CALL dsmsg( 5 )
 20   DO  i =1, maxpri
   ibasts = fcb( 2, i )
   IF ( ibasts == 0 ) CYCLE
   ibashi = ibasts + isysbf - 2
   ibaslo = ibasts - isysbf
   IF( indbas <= ibaslo .OR. indbas > ibashi ) CYCLE
   CALL dsmsg( 3 )
 END DO
 ibase( indbas )  = namfil
 fcb( 2, ifilex ) = indbas
 fcb(12, ifilex ) = indbas
 CALL dbmnam ( NAME, dname, ifilex )
 IF( iocode <= 1 ) GO TO 40
 IF( fcb( 13, ifilex ) == dname( 1 ) .AND.  &
     fcb( 14, ifilex ) == dname( 2 ) ) GO TO 35
!        CALL DBMSRF( DNAME, IUNI )
!        IF ( IUNI .EQ. IFILEX ) GO TO 35
 itrl(1) = NAME
 CALL rdtrl( itrl )
 DO  i = 2, 7
   IF ( itrl(i) /= 0 ) GO TO 35
 END DO
 IF ( iocode == 3 ) iocode = 1
 IF ( iocode == 2 ) iocode = 0
 GO TO 40
 35   CONTINUE
 nblock = fcb( 4,ifilex )
 IF ( nblock == 0 ) GO TO 40
 CALL dbmmgr ( 1 )
 indclr = fcb( 3, ifilex ) + indbas - 1
 indcbp = indclr
 GO TO 60
 40   nblock = 1
 fcb( 13, ifilex ) = dname( 1 )
 fcb( 14, ifilex ) = dname( 2 )
 CALL dbmmgr ( 1 )
 indclr = indbas + 5
 indcbp = indclr
 IF( iocode == 0 ) GO TO 60
 ibase( indbas+3 ) = 1
 ibase( indbas+4 ) = 0
 fcb( 8, ifilex )  = 0
 60   IF ( nblock == ibase( indbas+3 ) ) GO TO 70
 CALL dsmsg ( 102 )
 70   CALL dssdcb
!        PRINT *,' QOPEN,UN,CLR,BLK,IOP=',IFILEX,FCB(3,IFILEX),
!     &     FCB(4,IFILEX),IOP
 RETURN
END SUBROUTINE qopen
