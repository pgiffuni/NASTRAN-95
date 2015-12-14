SUBROUTINE dsxfsz
     INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 COMMON / xfist / lfist, nfist, ifist( 100 )
 COMMON / xfiat / ifiat( 643 )
 COMMON / zzzzzz/ mem(2)
 idsn   = ifilex
 nun    = 0
 itotal = 0
 10    lasblk = fcb( 6,idsn )
 ifrblk = fcb( 5,idsn )
 numblk = lasblk - ifrblk + 1
 IF ( idsn == ifilex ) GO TO 20
 nun    = nun + 1
 itotal = itotal + numblk
 GO TO 40
 20    ipblks = numblk
 IF ( fcb( 10, ifilex ) == 0 ) GO TO 40
 INDEX  = fcb( 10, ifilex )
 lblock = mem( INDEX+3 )
 ipblks = ipblks + lblock
 40    idsn   = IAND( mdsfcb( 3,idsn ), maskh2 )
 IF ( idsn /= 0 ) GO TO 10
 lim    = 2 * nfist
 DO  i = 1,lim,2
   IF ( NAME /= ifist( i ) ) CYCLE
   IF ( ifist( i+1 ) <= 0 ) EXIT
   indx   = ifist( i+1 )
   ifiat( indx+7 ) = ipblks * 2**16 + nun * 2**8
   ifiat( indx+8 ) = itotal * 2**16
   EXIT
 END DO
 70    CONTINUE
 maxusm = 0
 maxusd = 0
! ACCUMULATE TOTAL I/O USAGE STATISTICS
 DO  i = 1, 80
   IF ( i == 7 ) CYCLE
   itotl1 = 0
   itotl2 = 0
   IF ( fcb( 4, i ) == 0 ) CYCLE
   nexblk  = fcb( 10, i )
   IF ( nexblk /= 0 ) itotl1 = mem( nexblk+3 )
   IF ( fcb( 5, i ) /= 0 ) itotl2 = fcb( 6, i ) - fcb( 5, i ) + 1
   maxusm = maxusm + itotl1
   maxusd = maxusd + itotl2
 END DO
 IF ( maxblk < maxusm ) maxblk = maxusm
 IF ( maxdsk < maxusd ) maxdsk = maxusd
 RETURN
END SUBROUTINE dsxfsz
