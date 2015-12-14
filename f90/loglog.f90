SUBROUTINE loglog (a,b,c,d,e,f)
     
!     WRITTEN BY G.CHAN/UNISYS  7/92, THE 1992 SUMMER OLYMPIC WEEK
 
!     LOG-LOG TABLE LOOKUP           10 +------+------+------+--+
!                                     D                 *
!     INPUT : A,B, C,D, AND E         8 +------+-------/-----+--+
!     OUTPUT: F                                       /
!                                                    /
!     ALL A,B,C,D,E,F IN LOG          4 +------+----/-+------+--+
!     SCALE                           F            *
!                                                 /
!     LINEAR EVALUATION ON LOG        2 +------+-/-----+------+--+
!     SCALE (NO POLYNOMIAL            B         *
!     EVALUATION)
!                                     1 +------+------+------+--+
!                                       1      2A  E  4 C    8  10
 aa = ALOG10(a)
 bb = ALOG10(b)
 cc = ALOG10(c) - aa
 dd = ALOG10(d) - bb
 ee = ALOG10(e) - aa
 f  = 10.**(ee*dd/cc + bb)
 RETURN
 
 
 ENTRY smilog (a,b,c,d,e,f)
!     ==========================
 
!     SEMI-LOG TABLE LOOKUP       10 +--+--+--+--+--+--+--+--+--+
!                                  D                *
!     INPUT : A,B, C,D, AND E      8 +--+--+--+--+-/+--+--+--+--+
!     OUTPUT: F                                   /
!                                                /
!     A,C,E IN LINEAR SCALE        4 +--+--+--+-/+--+--+--+--+--+
!     B,D,F IN LOG SCALE           F           *
!                                             /
!                                  2 +--+--+-/+--+--+--+--+--+--+
!                                  B        *
 
!                                  1 +--+--+--+--+--+--+--+--+--+
!                                    0  1  2A 3E 4  C  6  7  8  9
 bb = ALOG10(b)
 cc = c - a
 dd = ALOG10(d) - bb
 ee = e - a
 f  = 10.**(ee*dd/cc + bb)
 RETURN
 
 
 ENTRY logsmi (a,b,c,d,e,f)
!     ==========================
 
!     LOG-SEMI TABLE LOOKUP          10 +-----+-----+-----+--+
!                                     D                 *
!     INPUT:  A,B, C,D, AND E         8 +-----+-----+--/--+--+
!     OUTPUT: F                                       /
!                                     6 +-----+-----+/----+--+
!     A,C,E IN LOG SCALE                            /
!     B,D,F IN LINEAR SCALE           4 +-----+----/+-----+--+
!                                     F           *
!                                     2 +-----+--/--+-----+--+
!                                     B         *
!                                     0 +-----+-----+-----+--+
!                                       1     2 A E 4   C 8  10
 aa = ALOG10(a)
 cc = ALOG10(c) - aa
 dd = d - b
 ee = ALOG10(e) - aa
 f  = ee*dd/cc + b
 RETURN
END SUBROUTINE loglog
