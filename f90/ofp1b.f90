SUBROUTINE ofp1b (line)
     
!     THIS SUBROUTINE WAS FORMED ONLY TO REDUCE THE SIZE OF OFP1 FOR
!     COMPILATION PURPOSES.  IT IS CALLED ONLY BY OFP1.
!     PREVIOUSLY THIS ROUTINE WAS NAMED OPF1A.
 
 
 INTEGER, INTENT(IN)                      :: line
 DIMENSION       fd(50),id(50),of(6)
 COMMON /system/ ibuf,l,nogo
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (core(1),of(1)),(id(1),fd(1),of(6))
 DATA            idum1, idum2, idum3, idum4, idum5, idum6 /  &
     4HDUM1,4HDUM2,4HDUM3,4HDUM4,4HDUM5,4HDUM6/, idum7, idum8, idum9 /  &
     4HDUM7,4HDUM8,4HDUM9/
 
 IF (line > 294) GO TO 10
 local = line - 174
 GO TO (175,176,177,178,179,180,181,182,183,184,185,186,187,188,189  &
     ,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204  &
     ,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219  &
     ,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234  &
     ,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249  &
     ,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264  &
     ,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279  &
     ,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294 ), local
 10 IF (line > 380) RETURN
 local = line - 294
 GO TO (295,296,297,298,299,300,301,302,303,304,305,306,307,308,309  &
     ,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324  &
     ,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339  &
     ,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354  &
     ,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369  &
     ,370,371,372,373,374,375,376,377,378,379,380) local
 175 CONTINUE
 GO TO 1000
 176 WRITE (l,676)
 GO TO 1000
 177 WRITE (l,677)
 GO TO 1000
 178 WRITE (l,678)
 GO TO 1000
 179 WRITE (l,679)
 GO TO 1000
 180 WRITE (l,680)
 GO TO 1000
 181 WRITE (l,681)
 GO TO 1000
 182 WRITE (l,682)
 GO TO 1000
 183 WRITE (l,683)
 GO TO 1000
 184 WRITE (l,684)
 GO TO 1000
 185 WRITE (l,685)
 GO TO 1000
 186 WRITE (l,686)
 GO TO 1000
 187 WRITE (l,687)
 GO TO 1000
 188 WRITE (l,688)
 GO TO 1000
 189 WRITE (l,689)
 GO TO 1000
 190 WRITE (l,690)
 GO TO 1000
 191 WRITE (l,691)
 GO TO 1000
 192 WRITE (l,692)
 GO TO 1000
 193 WRITE (l,693)
 GO TO 1000
 194 WRITE (l,694)
 GO TO 1000
 195 WRITE (l,695)
 GO TO 1000
 196 WRITE (l,696)
 GO TO 1000
 197 WRITE (l,697)
 GO TO 1000
 198 WRITE (l,698)
 GO TO 1000
 199 WRITE (l,699)
 GO TO 1000
 200 WRITE (l,700)
 GO TO 1000
 201 WRITE (l,701)
 GO TO 1000
 202 WRITE (l,702) id(3)
 GO TO 1000
 203 WRITE (l,703)
 GO TO 1000
 204 WRITE (l,704)
 GO TO 1000
 205 WRITE (l,705)
 GO TO 1000
 206 WRITE (l,706)
 GO TO 1000
 207 WRITE (l,707)
 GO TO 1000
 208 WRITE (l,708)
 GO TO 1000
 209 WRITE (l,709)
 GO TO 1000
 210 WRITE (l,710)
 GO TO 1000
 211 WRITE (l,711)
 GO TO 1000
 212 WRITE (l,712)
 GO TO 1000
 213 WRITE (l,713)
 GO TO 1000
 214 WRITE (l,714) id(5)
 GO TO 1000
 215 WRITE (l,715)
 GO TO 1000
 216 WRITE (l,716) (id(k),k=11,14),id(17),fd(18),(id(k),k=19,21)
 IF (id(17) == 1) WRITE (l,7161)
 IF (id(17) /= 1) WRITE (l,7162)
 IF (id(17) /= 1) nogo = 14
 GO TO 1000
 217 WRITE (l,717)
 GO TO 1000
 218 WRITE (l,718)
 GO TO 1000
 219 WRITE (l,719)
 GO TO 1000
 220 WRITE (l,720)
 GO TO 1000
 221 WRITE (l,721)
 GO TO 1000
 222 WRITE (l,722)
 GO TO 1000
 223 WRITE (l,723)
 GO TO 1000
 224 WRITE (l,724)
 GO TO 1000
 225 WRITE (l,725)
 GO TO 1000
 226 WRITE (l,726)
 GO TO 1000
 227 WRITE (l,727)
 GO TO 1000
 228 WRITE (l,728)
 GO TO 1000
 229 idd   = MOD(id(5),500000)
 jharm = (id(5)-idd)/500000
 iharm = (jharm-1)/2
 IF (MOD(jharm,2) == 1) GO TO 1229
 WRITE (l,729) idd,iharm
 GO TO 1000
 1229 WRITE (l,1729) idd,iharm
 GO TO 1000
 230 WRITE (l,730)
 GO TO 1000
 231 WRITE (l,731)
 GO TO 1000
 232 WRITE (l,732)
 GO TO 1000
 233 WRITE (l,733)
 GO TO 1000
 234 WRITE (l,734)
 GO TO 1000
 235 WRITE (l,735)
 GO TO 1000
 236 WRITE (l,736)
 GO TO 1000
 237 WRITE (l,737)
 GO TO 1000
 238 WRITE (l,738)
 GO TO 1000
 239 WRITE (l,739)
 GO TO 1000
 240 WRITE (l,740)
 GO TO 1000
 241 WRITE (l,741)
 GO TO 1000
 242 WRITE (l,742)
 GO TO 1000
 243 WRITE (l,743)
 GO TO 1000
 244 WRITE (l,744)
 GO TO 1000
 245 WRITE (l,745)
 GO TO 1000
 246 WRITE (l,746)
 GO TO 1000
 247 WRITE (l,747)
 GO TO 1000
 248 WRITE (l,748)
 GO TO 1000
 249 WRITE (l,749)
 GO TO 1000
 250 WRITE (l,750)
 GO TO 1000
 251 WRITE (l,751)
 GO TO 1000
 252 WRITE (l,752)
 GO TO 1000
 253 WRITE (l,753)
 GO TO 1000
 254 idx = idum1
 GO TO 1754
 255 idx = idum2
 GO TO 1754
 256 idx = idum3
 GO TO 1754
 257 idx = idum4
 GO TO 1754
 258 idx = idum5
 GO TO 1754
 259 idx = idum1
 GO TO 1759
 260 idx = idum2
 GO TO 1759
 261 idx = idum3
 GO TO 1759
 262 idx = idum4
 GO TO 1759
 263 idx = idum5
 GO TO 1759
 264 WRITE (l,764)
 GO TO 1000
 265 WRITE (l,765)
 GO TO 1000
 266 idx = idum1
 GO TO 1766
 267 idx = idum2
 GO TO 1766
 268 idx = idum3
 GO TO 1766
 269 idx = idum4
 GO TO 1766
 270 idx = idum5
 GO TO 1766
 271 idx = idum1
 GO TO 1771
 272 idx = idum2
 GO TO 1771
 273 idx = idum3
 GO TO 1771
 274 idx = idum4
 GO TO 1771
 275 idx = idum5
 GO TO 1771
 276 WRITE (l,776)
 GO TO 1000
 277 WRITE (l,777)
 GO TO 1000
 278 WRITE (l,778)
 GO TO 1000
 279 WRITE (l,779)
 GO TO 1000
 280 idx = idum6
 GO TO 1754
 281 idx = idum7
 GO TO 1754
 282 idx = idum8
 GO TO 1754
 283 idx = idum9
 GO TO 1754
 284 idx = idum6
 GO TO 1759
 285 idx = idum7
 GO TO 1759
 286 idx = idum8
 GO TO 1759
 287 idx = idum9
 GO TO 1759
 288 idx = idum6
 GO TO 1766
 289 idx = idum7
 GO TO 1766
 290 idx = idum8
 GO TO 1766
 291 idx = idum9
 GO TO 1766
 292 idx = idum6
 GO TO 1771
 293 idx = idum7
 GO TO 1771
 294 idx = idum8
 GO TO 1771
 295 idx = idum9
 GO TO 1771
 296 WRITE (l,796)
 GO TO 1000
 297 WRITE (l,797)
 GO TO 1000
 298 WRITE (l,798)
 GO TO 1000
 299 WRITE (l,799)
 GO TO 1000
 300 WRITE (l,800)
 GO TO 1000
 301 WRITE (l,801)
 GO TO 1000
 302 WRITE (l,802)
 GO TO 1000
 303 WRITE (l,803)
 GO TO 1000
 304 WRITE (l,804)
 GO TO 1000
 305 WRITE (l,805)
 GO TO 1000
 306 WRITE (l,806)
 GO TO 1000
 307 WRITE (l,807)
 GO TO 1000
 308 WRITE (l,808)
 GO TO 1000
 309 WRITE (l,809)
 GO TO 1000
 310 WRITE (l,810)
 GO TO 1000
 311 WRITE (l,811)
 GO TO 1000
 312 WRITE (l,812)
 GO TO 1000
 313 WRITE (l,813)
 GO TO 1000
 314 WRITE (l,814)
 GO TO 1000
 315 WRITE (l,815)
 GO TO 1000
 316 WRITE (l,816)
 GO TO 1000
 317 WRITE (l,817)
 GO TO 1000
 318 WRITE (l,818)
 GO TO 1000
 319 WRITE (l,819)
 GO TO 1000
 320 WRITE (l,820)
 GO TO 1000
 321 WRITE (l,821)
 GO TO 1000
 322 WRITE (l,822)
 GO TO 1000
 323 WRITE (l,823)
 GO TO 1000
 324 WRITE (l,824)
 GO TO 1000
 325 WRITE (l,825)
 GO TO 1000
 326 WRITE (l,826)
 GO TO 1000
 327 WRITE (l,827)
 GO TO 1000
 328 idd = id(3) - 64
 WRITE (l,828) idd
 GO TO 1000
 329 WRITE (l,829)
 GO TO 1000
 330 WRITE (l,830)
 GO TO 1000
 331 idd = id(3) - 64
 WRITE (l,831) idd
 GO TO 1000
 332 WRITE (l,832)
 GO TO 1000
 333 WRITE (l,833)
 GO TO 1000
 334 WRITE (l,834)
 GO TO 1000
 335 WRITE (l,835)
 GO TO 1000
 336 WRITE (l,836)
 GO TO 1000
 337 WRITE (l,837)
 GO TO 1000
 338 WRITE (l,838)
 GO TO 1000
 339 WRITE (l,839)
 GO TO 1000
 340 WRITE (l,840)
 GO TO 1000
 341 WRITE (l,841)
 GO TO 1000
 342 WRITE (l,842)
 GO TO 1000
 343 WRITE (l,843)
 GO TO 1000
 344 WRITE (l,844)
 GO TO 1000
 345 WRITE (l,845)
 GO TO 1000
 346 WRITE (l,846)
 GO TO 1000
 347 WRITE (l,847)
 GO TO 1000
 348 WRITE (l,848)
 GO TO 1000
 349 WRITE (l,849)
 GO TO 1000
 350 WRITE (l,850)
 GO TO 1000
 351 WRITE (l,851)
 GO TO 1000
 352 WRITE (l,852)
 GO TO 1000
 353 WRITE (l,853)
 GO TO 1000
 354 WRITE (l,854) id(6),id(7),fd(3)
 GO TO 1000
 355 WRITE (l,855)
 GO TO 1000
 356 WRITE (l,856)
 GO TO 1000
 357 WRITE (l,857)
 GO TO 1000
 358 WRITE (l,858)
 GO TO 1000
 359 WRITE (l,859)
 GO TO 1000
 360 WRITE (l,860)
 GO TO 1000
 361 WRITE (l,861)
 GO TO 1000
 362 WRITE (l,862)
 GO TO 1000
 363 WRITE (l,863)
 GO TO 1000
 364 WRITE (l,864)
 GO TO 1000
 365 WRITE (l,865)
 GO TO 1000
 366 WRITE (l,866)
 GO TO 1000
 367 WRITE (l,867)
 GO TO 1000
 368 WRITE (l,868)
 GO TO 1000
 369 WRITE (l,869)
 GO TO 1000
 370 WRITE (l,870)
 GO TO 1000
 371 WRITE (l,871)
 GO TO 1000
 372 WRITE (l,872)
 GO TO 1000
 373 WRITE (l,873)
 GO TO 1000
 374 WRITE (l,874)
 GO TO 1000
 375 WRITE (l,875)
 GO TO 990
 376 WRITE (l,876)
 GO TO 1000
 377 WRITE (l,877)
 GO TO 1000
 378 WRITE (l,878)
 GO TO 1000
 379 WRITE (l,879)
 GO TO 1000
 380 WRITE (l,880)
 GO TO 1000
 990 WRITE  (l,995)
 995 FORMAT (1H )
 1000 RETURN
 
!     ******************************************************************
 
  676 FORMAT (2(25X,5HAXIAL,30X), /2(7X,4HTIME,14X,5HFORCE,9X,6HTORQUE, &
            15X))
  677 FORMAT (21X,17HBEND-MOMENT-END-A,12X,17HBEND-MOMENT-END-B,18X, &
             5HSHEAR,17X, /7X,4HTIME,3(8X,7HPLANE 1,7X,7HPLANE 2),9X, &
             5HFORCE,10X,6HTORQUE)
  678 FORMAT (2(25X,5HFORCE,10X,5HFORCE,15X), /2(7X,4HTIME,13X, &
             7HPTS 1,3,8X,7HPTS 2,4,14X))
  679 FORMAT (2(24X,6HMOMENT,9X,6HMOMENT,15X), /2(7X,4HTIME,13X, &
             7HPTS 1,3,8X,7HPTS 2,4,14X))
  680 FORMAT (8X,4HTIME,3X,2(11X,11HBEND-MOMENT),11X,12HTWIST-MOMENT, &
             13X,5HSHEAR,17X,5HSHEAR, /31X,1HX,21X,1HY,43X,1HX,21X,1HY)
  681 FORMAT (4(8X,4HTIME,10X,5HFORCE,6X))
  682 FORMAT (2(21X,5HAXIAL,7X,6HSAFETY,6X,9HTORSIONAL,5X,6HSAFETY), / &
             2(7X,4HTIME,9X,6HSTRESS,7X,6HMARGIN,8X,6HSTRESS,6X, &
             6HMARGIN))
  683 FORMAT (7X,4HTIME,12X,3HSA1,12X,3HSA2,12X,3HSA3,10X, &
             12HAXIAL-STRESS,8X,6HSA-MAX,9X,6HSA-MIN,11X,6HM.S.-T, /23X, &
             3HSB1,12X,3HSB2,12X,3HSB3,30X,6HSB-MAX,9X,6HSB-MIN,11X, &
             6HM.S.-C)
  684 FORMAT (2(26X,7HMAXIMUM,8X,7HAVERAGE,6X,6HSAFETY), /2(8X,4HTIME, &
             15X,5HSHEAR,10X,5HSHEAR,7X,6HMARGIN))
  685 FORMAT (2(54X,6HSAFETY), /2(7X,4HTIME,15X,7HMAXIMUM,8X,7HAVERAGE, &
             6X,6HMARGIN))
  686 FORMAT (19X,5HFIBRE,11X,32HSTRESSES IN ELEMENT COORD SYSTEM,13X, &
             31HPRINCIPAL STRESSES (ZERO SHEAR),10X,7HMAXIMUM, /7X, &
             4HTIME,7X,8HDISTANCE,7X,8HNORMAL-X,7X,8HNORMAL-Y,6X, &
             8HSHEAR-XY,7X,5HANGLE,9X,5HMAJOR,11X,5HMINOR,10X,5HSHEAR)
  687 FORMAT (20X,32HSTRESSES IN ELEMENT COORD SYSTEM,12X,9HPRINCIPAL, &
             11X,18HPRINCIPAL STRESSES,10X,7HMAXIMUM, /7X,4HTIME,8X, &
             8HNORMAL-X,6X,8HNORMAL-Y,7X,8HSHEAR-XY,6X,12HSTRESS ANGLE, &
             9X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
  688 FORMAT (4(8X,4HTIME, 9X,6HSTRESS,6X))
  689 FORMAT (5X,4HTIME,15X,3HSA1,12X,3HSA2,12X,3HSA3,12X,3HSA4,8X, &
             12HAXIAL-STRESS,6X,6HSA-MAX,9X,6HSA-MIN,5X,6HM.S.-T, /24X, &
             3HSB1,12X,3HSB2,12X,3HSB3,12X,3HSB4,26X,6HSB-MAX,9X, &
             6HSB-MIN,5X,6HM.S.-C)
  690 FORMAT (53X,5HAXIAL, /13X,9HFREQUENCY,31X,5HFORCE,41X,6HTORQUE)
  691 FORMAT (11X,2(42X,5HFORCE), /13X,9HFREQUENCY,30X,7HPTS 1,3,40X, &
             7HPTS 2,4)
  692 FORMAT (11X,2(41X,6HMOMENT), /13X,9HFREQUENCY,30X,7HPTS 1,3,40X, &
             7HPTS 2,4)
  693 FORMAT (5X,9HFREQUENCY,2X,2(11X,11HBEND-MOMENT),10X, &
             12HTWIST-MOMENT,2(13X,5HSHEAR,4X), /2(31X,1HX,21X,1HY,12X))
  694 FORMAT (2(12X,9HFREQUENCY,20X,5HFORCE,12X))
  695 FORMAT (53X,5HAXIAL,39X,9HTORSIONAL, /13X,9HFREQUENCY, &
             2(30X,6HSTRESS,11X))
  696 FORMAT (52X,7HMAXIMUM,39X,7HAVERAGE, /13X,9HFREQUENCY, &
             2(31X,5HSHEAR,10X))
  697 FORMAT (20X,5HFIBRE,37X,'- STRESSES IN ELEMENT COORDINATE SYSTEM', &
             2H -, /4X,'FREQUENCY,6X,8HDISTANCE',18X,8HNORMAL-X,26X, &
             8HNORMAL-Y,25X,8HSHEAR-XY)
  698 FORMAT (53X,41H- STRESSES IN ELEMENT COORDINATE SYSTEM -, /9X, &
             9HFREQUENCY,18X,8HNORMAL-X,26X,8HNORMAL-Y,26X,8HSHEAR-XY)
  699 FORMAT (2(12X,9HFREQUENCY,19X,6HSTRESS,12X))
  700 FORMAT (39X,4(8HLOCATION,7X),6X,7HAVERAGE, /8X,9HFREQUENCY,26X, &
             1H1,14X,1H2,14X,1H3,14X,1H4,13X,12HAXIAL STRESS)
  701 FORMAT (21X,17HBEND-MOMENT-END-A,12X,17HBEND-MOMENT-END-B,18X, &
             5HSHEAR,17X, /4X,9HFREQUENCY,3(6X,7HPLANE 1,7X, &
             9HPLANE 2  ),6X,5HFORCE,10X,6HTORQUE)
  702 FORMAT (27X,'O U T P U T   F R O M   G R I D   P O I N T   W E I', &
             ' G H T   G E N E R A T O R', /1H0,53X, &
             17HREFERENCE POINT =,I9)
  703 FORMAT (5X,9HSECTOR-ID,/6X,8HPOINT-ID,/7X,7HRING-ID,2X,8HHARMONIC, &
             8X,2HT1,13X,2HT2,13X,2HT3,13X,2HR1,13X,2HR2,13X,2HR3)
  704 FORMAT (11X,'S T R E S S E S   I N   A X I S - S Y M M E T R I C', &
             '   C O N I C A L   S H E L L   E L E M E N T S   ', &
             '(CCONEAX)')
  705 FORMAT (13X,'F O R C E S   I N   A X I S - S Y M M E T R I C   ', &
             'C O N I C A L   S H E L L   E L E M E N T S   (CCONEAX)')
  706 FORMAT (8H ELEMENT,10X,5HPOINT,5X,5HFIBRE,11X,'STRESSES IN ELEM', &
             'ENT COORD SYSTEM',8X,'PRINCIPAL STRESSES (ZERO SHEAR)', &
             8X,7HMAXIMUM, /3X,'ID.  HARMONIC  ANGLE    DISTANCE',7X, &
             8HNORMAL-V,6X,8HNORMAL-U,6X,8HSHEAR-UV,6X,5HANGLE,7X, &
             5HMAJOR,9X,5HMINOR,9X,5HSHEAR)
  707 FORMAT (9H  ELEMENT,5X,8HHARMONIC,4X,5HPOINT,4X,2(7X,5HBEND-, &
             6HMOMENT),6X,12HTWIST-MOMENT,2(11X,5HSHEAR,1X), /3X,3HID., &
             9X,6HNUMBER,5X,5HANGLE,15X,1HV,17X,1HU,37X,1HV,16X,1HU)
  708 FORMAT (31X,'C O M P L E X   D I S P L A C E M E N T   ', &
             'V E C T O R  (SOLUTION SET)')
  709 FORMAT (35X,'C O M P L E X   V E L O C I T Y   V E C T O R  ', &
             '(SOLUTION SET)')
  710 FORMAT (31X,'C O M P L E X   A C C E L E R A T I O N   ', &
             'V E C T O R   (SOLUTION SET)')
  711 FORMAT (43X,46HV E L O C I T Y   V E C T O R   (SOLUTION SET)) 
  712 FORMAT (39X,'D I S P L A C E M E N T   V E C T O R   ', &
             '(SOLUTION SET)')
  713 FORMAT (39X,'A C C E L E R A T I O N   V E C T O R   ', &
             '(SOLUTION SET)')
  714 FORMAT (29X,'C O M P L E X   E I G E N V E C T O R   N O .',I11, &
             3X,14H(SOLUTION SET))
  715 FORMAT (30X,'E I G E N V A L U E   A N A L Y S I S   ', &
             'S U M M A R Y   (GIVENS METHOD)')
  716 FORMAT (///,36X,45HNUMBER OF EIGENVALUES EXTRACTED . . . . . . ., &
             I10,//36X,45HNUMBER OF EIGENVECTORS COMPUTED . . . . . . .,&
             I10,//36X,45HNUMBER OF EIGENVALUE CONVERGENCE FAILURES . .,&
             I10,//36X,45HNUMBER OF EIGENVECTOR CONVERGENCE FAILURES. .,&
            I10,///36X,45HREASON FOR TERMINATION. . . . . . . . . . . .,&
        I10,1H*,///36X,45HLARGEST OFF-DIAGONAL MODAL MASS TERM. . . . .,&
       1P,E10.2,//76X,5H. . .,I10,/46X,'MODE PAIR. . . . . . . . . . .',&
           /76X,5H. . .,I10,//36X,33HNUMBER OF OFF-DIAG0NAL MODAL MASS ,&
           /41X,40HTERMS FAILING CRITERION. . . . . . . . .,I10)
 7161 FORMAT (//36X,22H(* NORMAL TERMINATION))
 7162 FORMAT (//36X,31H(* INSUFFICIENT TIME REMAINING))
  717 FORMAT (107X,22HOCTAHEDRAL    PRESSURE, /6X,10HELEMENT-ID,8X, &
             8HSIGMA-XX,6X,8HSIGMA-YY,6X,8HSIGMA-ZZ,7X,6HTAU-YZ,8X, &
             6HTAU-XZ,8X,6HTAU-XY,8X,5HTAU-0,10X,1HP)
  718 FORMAT (107X,22HOCTAHEDRAL    PRESSURE, /6X,10H TIME     ,8X, &
             8HSIGMA-XX,6X,8HSIGMA-YY,6X,8HSIGMA-ZZ,7X,6HTAU-YZ,8X, &
             6HTAU-XZ,8X,6HTAU-XY,8X,5HTAU-0,10X,1HP)
  719 FORMAT (18X,10HELEMENT-ID,8X,8HSIGMA-XX,6X,8HSIGMA-YY,6X, &
             8HSIGMA-ZZ,7X,6HTAU-YZ,8X,6HTAU-XZ,8X,6HTAU-XY)
  720 FORMAT (18X,10HFREQUENCY ,8X,8HSIGMA-XX,6X,8HSIGMA-YY,6X, &
             8HSIGMA-ZZ,7X,6HTAU-YZ,8X,6HTAU-XZ,8X,6HTAU-XY)
  721 FORMAT (19X,'S T R E S S E S   I N   S O L I D   T E T R A H E D', &
             ' R O N   E L E M E N T S   ( C T E T R A )')
  722 FORMAT (11X,'C O M P L E X   S T R E S S E S   I N   S O L I D  ', &
           ' T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )')
  723 FORMAT (25X,'S T R E S S E S   I N   S O L I D   W E D G E   ', &
             'E L E M E N T S   ( C W E D G E )')
  724 FORMAT (17X,'C O M P L E X   S T R E S S E S   I N   S O L I D  ', &
             ' W E D G E   E L E M E N T S   ( C W E D G E )')
  725 FORMAT (20X,'S T R E S S E S   I N   S O L I D   H E X A H E D R', &
             ' O N   E L E M E N T S   ( C H E X A 1 )')
  726 FORMAT (12X,'C O M P L E X   S T R E S S E S   I N   S O L I D  ', &
             ' H E X A H E D R O N   E L E M E N T S   ( C H E X A 1 )')
  727 FORMAT (20X,'S T R E S S E S   I N   S O L I D   H E X A H E D R', &
             ' O N   E L E M E N T S   ( C H E X A 2 )')
  728 FORMAT (12X,'C O M P L E X   S T R E S S E S   I N   S O L I D  ', &
             ' H E X A H E D R O N   E L E M E N T S   ( C H E X A 2 )')
  729 FORMAT (6X,10HPOINT-ID =,I7,4X,10HHARMONIC =,I4) 
 1729 FORMAT (6X,10HPOINT-ID =,I7,4X,10HHARMONIC =,I4,1H*)
  730 FORMAT (5X,8HHARMONIC,5(3X,8HPOINT-ID,5X,2HT1,5X))
  731 FORMAT (10X,'V E L O C I T I E S   I N   A X I S Y M M E T R I C', &
             '   F L U I D   E L E M E N T S   ( C A X I F 2 - ', &
             'S T R E S S )')
  732 FORMAT (10X,'V E L O C I T I E S   I N   A X I S Y M M E T R I C', &
             '   F L U I D   E L E M E N T S   ( C A X I F 3 - ', &
             'S T R E S S )')
  733 FORMAT (10X,'V E L O C I T I E S   I N   A X I S Y M M E T R I C', &
             '   F L U I D   E L E M E N T S   ( C A X I F 4 - ', &
             'S T R E S S )')
  734 FORMAT (24X,'V E L O C I T I E S   I N   S L O T   E L E M E N T', &
             ' S   ( C S L O T 3 - S T R E S S )')
  735 FORMAT (24X,'V E L O C I T I E S   I N   S L O T   E L E M E N T', &
             ' S   ( C S L O T 4 - S T R E S S )')
  736 FORMAT (2X,'C O M P L E X   V E L O C I T I E S   I N   A X I S ', &
             'Y M M E T R I C   F L U I D   E L E M E N T S   ', &
             '( C A X I F 2 - S T R E S S )')
  737 FORMAT (2X,'C O M P L E X   V E L O C I T I E S   I N   A X I S ', &
             'Y M M E T R I C   F L U I D   E L E M E N T S   ', &
             '( C A X I F 3 - S T R E S S )')
  738 FORMAT (2X,'C O M P L E X   V E L O C I T I E S   I N   A X I S ', &
             'Y M M E T R I C   F L U I D   E L E M E N T S   ', &
             '( C A X I F 4 - S T R E S S )')
  739 FORMAT (15X,'C O M P L E X   V E L O C I T I E S   I N   S L O T', &
             '  E L E M E N T S   ( C S L O T 3 - S T R E S S )')
  740 FORMAT (15X,'C O M P L E X   V E L O C I T I E S   I N   S L O T',&
             '  E L E M E N T S   ( C S L O T 4 - S T R E S S )')
  741 FORMAT (8X,7HELEMENT,17X,6HCENTER,25X,7HEDGE  1,19X,7HEDGE  2,19X,&
             7HEDGE  3,/10X,3HID., 8X,27HR --------- PHI --------- Z,&
             12X,15HS --------- PHI,11X,15HS --------- PHI,11X,&
             15HS --------- PHI)
  742 FORMAT (31X,6HCENTER,25X,7HEDGE  1,19X,7HEDGE  2,19X,7HEDGE  3,&
             /2X,8H   TIME ,10X,27HR --------- PHI --------- Z,12X,&
             15HS --------- PHI,11X,15HS --------- PHI,11X,&
             15HS --------- PHI)
  743 FORMAT (32X,6HCENTER,25X,7HEDGE  1,19X,7HEDGE  2,19X,7HEDGE  3,&
             /4X,9HFREQUENCY,8X,27HR --------- PHI --------- Z,12X,&
             15HS --------- PHI,11X,15HS --------- PHI,11X,&
             15HS --------- PHI)
  744 FORMAT (13X,7HELEMENT,18X,6HCENTER,20X,7HEDGE  1,11X,7HEDGE  2,&
             11X,7HEDGE  3,11X,7HEDGE  4,/15X,3HID.,13X,&
             19HR --------------- Z,17X,1HS,17X,1HS,17X,1HS,17X,1HS)
  745 FORMAT (38X,6HCENTER,20X,7HEDGE  1,11X,7HEDGE  2,11X,7HEDGE  3,&
             11X,7HEDGE  4, /11X,4HTIME,16X,19HR --------------- Z,17X,&
             1HS,17X,1HS,17X,1HS,17X,1HS)
  746 FORMAT (38X,6HCENTER,20X,7HEDGE  1,11X,7HEDGE  2,11X,7HEDGE  3,&
             11X,7HEDGE  4,/9X,9HFREQUENCY,13X,19HR --------------- Z,&
             17X,1HS,17X,1HS,17X,1HS,17X,1HS)
  747 FORMAT (9X,7HELEMENT,24X,6HCENTER,26X,7HEDGE  1,15X,7HEDGE  2,&
             15X,7HEDGE  3,/11X,3HID.,17X,23HR ------------------- Z,&
             21X,1HS,21X,1HS,21X,1HS)
  748 FORMAT (40X,6HCENTER,26X,7HEDGE  1,15X,7HEDGE  2,15X,7HEDGE  3,&
             /7X,4HTIME,20X,23HR ------------------- Z,21X,1HS,21X,1HS,&
             21X,1HS)
  749 FORMAT (40X,6HCENTER,26X,7HEDGE  1,15X,7HEDGE  2,15X,7HEDGE  3,&
             /5X,9HFREQUENCY,17X,23HR ------------------- Z,21X,1HS,21X,&
             1HS,21X,1HS)
  750 FORMAT (14X,7HELEMENT,30X,6HCENTER,47X,4HEDGE,/16X,3HID.,21X,&
             27HR ----------------------- Z,25X,&
             28HS ---------------------- PHI)
  751 FORMAT (51X,6HCENTER,47X,4HEDGE, /12X,4HTIME,24X,&
             27HR ----------------------- Z,25X,&
             28HS ---------------------- PHI)
  752 FORMAT (51X,6HCENTER,47X,4HEDGE, /10X,9HFREQUENCY,21X,&
             27HR ----------------------- Z,25X,&
             28HS ---------------------- PHI)
  753 FORMAT (46X,35HT E M P E R A T U R E   V E C T O R)
 1754 WRITE  (L,754) IDX
      GO TO  1000
  754 FORMAT (36X,'S T R E S S E S   I N   U S E R   E L E M E N T S',&
             '  (C',A4,1H))
 1759 WRITE  (L,759) IDX
      GO TO  1000
  759 FORMAT (38X,'F O R C E S   I N   U S E R   E L E M E N T S   (C',&
             A4,1H) )
  764 FORMAT (5X,9H    EL-ID,6X,2HS1,11X,2HS2,11X,2HS3,11X,2HS4,11X,&
             2HS5,11X,2HS6,11X,2HS7,11X,2HS8,11X,2HS9)
  765 FORMAT (5X,9H    EL-ID,6X,2HF1,11X,2HF2,11X,2HF3,11X,2HF4,11X,&
             2HF5,11X,2HF6,11X,2HF7,11X,2HF8,11X,2HF9)
 1766 WRITE  (L,766) IDX
      GO TO  1000
  766 FORMAT (28X,'C O M P L E X   S T R E S S E S   I N   U S E R   ',&
             'E L E M E N T S   (C',A4,1H))
 1771 WRITE  (L,771) IDX
      GO TO  1000
  771 FORMAT (30X,'C O M P L E X   F O R C E S   I N   U S E R   E L E',&
             ' M E N T S   (C',A4,1H))
  776 FORMAT (5X,9H     TIME,6X,2HS1,11X,2HS2,11X,2HS3,11X,2HS4,11X,&
             2HS5,11X,2HS6,11X,2HS7,11X,2HS8,11X,2HS9)
  777 FORMAT (5X,9H     TIME,6X,2HF1,11X,2HF2,11X,2HF3,11X,2HF4,11X,&
             2HF5,11X,2HF6,11X,2HF7,11X,2HF8,11X,2HF9)
  778 FORMAT (5X,9HFREQUENCY,6X,2HS1,11X,2HS2,11X,2HS3,11X,2HS4,11X,&
             2HS5,11X,2HS6,11X,2HS7,11X,2HS8,11X,2HS9)
  779 FORMAT (5X,9HFREQUENCY,6X,2HF1,11X,2HF2,11X,2HF3,11X,2HF4,11X,&
             2HF5,11X,2HF6,11X,2HF7,11X,2HF8,11X,2HF9)
  796 FORMAT (6X,'POINT ID.   TYPE',6X,'ID   VALUE     ID+1 VALUE    ',&
             'ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE')
  797 FORMAT (19X,'F I N I T E   E L E M E N T   T E M P E R A T U R E',&
             '   G R A D I E N T S   A N D   F L U X E S')
  798 FORMAT (4X,'ELEMENT-ID   EL-TYPE        X-GRADIENT       Y-',&
             'GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX',&
             '           Z-FLUX')
  799 FORMAT (4X,'TIME         EL-TYPE        X-GRADIENT       Y-',&
             'GRADIENT       Z-GRADIENT        X-FLUX           Y-FLUX',&
             '           Z-FLUX')
  800 FORMAT (26X,'ELEMENT-ID      APPLIED-LOAD       CONVECTION      ',&
             ' RADIATION           TOTAL')
  801 FORMAT (26X,'TIME            APPLIED-LOAD       CONVECTION      ',&
             ' RADIATION           TOTAL')
  802 FORMAT (33X,'H E A T   F L O W   I N T O   H B D Y   E L E M E N',&
             ' T S   (CHBDY)')
  803 FORMAT (6X,16HTIME        TYPE  ,6X,7H  VALUE)
  804 FORMAT (21X,'S T R E S S E S   I N   Q U A D R I L A T E R A L',&
             '   M E M B R A N E S      ( C Q D M E M 1 )')
  805 FORMAT (14X,'C O M P L E X   S T R E S S E S   I N   Q U A D R I',&
             ' L A T E R A L   M E M B R A N E S   ( C Q D M E M 1 )')
  806 FORMAT (26X,'S T R E S S E S   A C T I N G   I N   Q D M E M 2  ',&
             ' E L E M E N T S   (CQDMEM2)')
  807 FORMAT (19X,'C O M P L E X   S T R E S S E S   A C T I N G   I N',&
             '   Q D M E M 2   E L E M E N T S   (CQDMEM2)')
  808 FORMAT (18X,'S T R E S S E S   I N   G E N E R A L   Q U A D R I',&
             ' L A T E R A L   E L E M E N T S', 6X,15H( C Q U A D 4 ))
  809 FORMAT ('0*** THIS FORMAT 809/OFP1B NOT USED ***')
!                   ==============================
  810 FORMAT (28X,'F O R C E S   A C T I N G   O N   Q D M E M 2   E L',&
             ' E M E N T S   (CQDMEM2)')
  811 FORMAT (20X,'C O M P L E X   F O R C E S   A C T I N G   O N   ',&
             'Q D M E M 2   E L E M E N T S   (CQDMEM2)')
  812 FORMAT (18X,106H====== POINT  1 ======      ====== POINT  2 ======&
           &      ====== POINT  3 ======      ====== POINT  4 ====== , /7X,&
             7HELEMENT,4X,8HF-FROM-4,6X,8HF-FROM-2,6X,8HF-FROM-1,6X,&
             8HF-FROM-3,6X,8HF-FROM-2,6X,8HF-FROM-4,6X,8HF-FROM-3,6X,&
             8HF-FROM-1, /9X,2HID,15X,6HKICK-1,7X,8HSHEAR-12,7X,&
             6HKICK-2,7X,8HSHEAR-23,7X,6HKICK-3,7X,8HSHEAR-34,7X,&
             6HKICK-4,7X,8HSHEAR-41 )
  813 FORMAT (18X,106H====== POINT  1 ======      ====== POINT  2 ======&
            &      ====== POINT  3 ======      ====== POINT  4 ====== , /14X,&
             4X,8HF-FROM-4,6X,8HF-FROM-2,6X,8HF-FROM-1,6X,8HF-FROM-3,6X,&
             8HF-FROM-2,6X,8HF-FROM-4,6X,8HF-FROM-3,6X,8HF-FROM-1, /5X,&
             9HFREQUENCY,12X,6HKICK-1,7X,8HSHEAR-12,7X,6HKICK-2,7X,&
             8HSHEAR-23,7X,6HKICK-3,7X,8HSHEAR-34,7X,6HKICK-4,7X,&
             8HSHEAR-41)
  814 FORMAT (18X,106H====== POINT  1 ======      ====== POINT  2 ======&
            &      ====== POINT  3 ======      ====== POINT  4 ====== , /14X,&
             4X,8HF-FROM-4,6X,8HF-FROM-2,6X,8HF-FROM-1,6X,8HF-FROM-3,6X,&
             8HF-FROM-2,6X,8HF-FROM-4,6X,8HF-FROM-3,6X,8HF-FROM-1, /10X, &
             4HTIME,12X,6HKICK-1,7X,8HSHEAR-12,7X,6HKICK-2,7X,8HSHEAR-23, &
             7X,6HKICK-3,7X,8HSHEAR-34,7X,6HKICK-4,7X,8HSHEAR-41)
  815 FORMAT (6X,16HSUBCASE     TYPE,10X,2HT1,13X,2HT2,13X,2HT3,13X,&
             2HR1,13X,2HR2,13X,2HR3)
  816 FORMAT (2(26X,7HMAXIMUM,8X,7HAVERAGE,6X,6HSAFETY),/2(6X,7HSUBCASE,&
             14X,5HSHEAR,10X,5HSHEAR,7X,6HMARGIN))
  817 FORMAT (20X,32HSTRESSES IN ELEMENT COORD SYSTEM,12X,9HPRINCIPAL,&
             11X,18HPRINCIPAL STRESSES,10X,7HMAXIMUM, /6X,7HSUBCASE,6X,&
             8HNORMAL-X,6X,8HNORMAL-Y,7X,8HSHEAR-XY,6X,12HSTRESS ANGLE,&
             9X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
  818 FORMAT (6X,7HSUBCASE,11X,3HSA1,12X,3HSA2,12X,3HSA3,12X,3HSA4,8X,&
             12HAXIAL-STRESS,6X,6HSA-MAX,9X,6HSA-MIN,5X,6HM.S.-T, /24X,&
             3HSB1,12X,3HSB2,12X,3HSB3,12X,3HSB4,26X,6HSB-MAX,9X,&
             6HSB-MIN,5X,6HM.S.-C)
  819 FORMAT (18X,106H====== POINT  1 ======      ====== POINT  2 ======&
            &      ====== POINT  3 ======      ====== POINT  4 ====== , /14X,&
             4X,8HF-FROM-4,6X,8HF-FROM-2,6X,8HF-FROM-1,6X,8HF-FROM-3,6X,&
             8HF-FROM-2,6X,8HF-FROM-4,6X,8HF-FROM-3,6X,8HF-FROM-1, /5X,&
             7HSUBCASE,14X,6HKICK-1,7X,8HSHEAR-12,7X,6HKICK-2,7X,&
             8HSHEAR-23,7X,6HKICK-3,7X,8HSHEAR-34,7X,6HKICK-4,7X,&
             8HSHEAR-41 )
  820 FORMAT (21X,17HBEND-MOMENT-END-A,12X,17HBEND-MOMENT-END-B,18X,&
             5HSHEAR, /6X,7HSUBCASE,6X,3(7HPLANE 1,7X,7HPLANE 2,8X),&
             6H FORCE,10X,6HTORQUE)
  821 FORMAT (2(21X,5HAXIAL,7X,6HSAFETY,6X,9HTORSIONAL,5X,6HSAFETY), /&
             2(6X,7HSUBCASE,7X,6HSTRESS,7X,6HMARGIN,8X,6HSTRESS,6X,&
             6HMARGIN))
  822 FORMAT (2(25X,5HAXIAL,30X) , /2(6X,7HSUBCASE,12X,5HFORCE, 9X,&
             6HTORQUE,15X))
  823 FORMAT (5X,7HELEMENT,8X,33HSTRESSES IN MATERIAL COORD SYSTEM,12X,&
             9HPRINCIPAL,11X,18HPRINCIPAL STRESSES,12X,3HMAX)
  824 FORMAT (13X,7HELEMENT,33X,'- STRESSES IN MATERIAL COORDINATE ',&
             'SYSTEM -', /15X,3HID.,18X,8HNORMAL-X,26X,8HNORMAL-Y,26X,&
             8HSHEAR-XY)
  825 FORMAT (20X,33HSTRESSES IN MATERIAL COORD SYSTEM,11X,9HPRINCIPAL,&
             11X,18HPRINCIPAL STRESSES,10X,7HMAXIMUM, /7X,4HTIME,8X,&
             8HNORMAL-X,6X,8HNORMAL-Y,7X,8HSHEAR-XY,6X,12HSTRESS ANGLE,&
             9X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
  826 FORMAT (53X,42H- STRESSES IN MATERIAL COORDINATE SYSTEM -, /9X,&
             9HFREQUENCY,18X,8HNORMAL-X,26X,8HNORMAL-Y,26X,8HSHEAR-XY)
  827 FORMAT (20X,33HSTRESSES IN MATERIAL COORD SYSTEM,11X,9HPRINCIPAL,&
             11X,18HPRINCIPAL STRESSES,10X,7HMAXIMUM, /6X,7HSUBCASE,6X,&
             8HNORMAL-X,6X,8HNORMAL-Y,7X,8HSHEAR-XY,6X,12HSTRESS ANGLE,&
             9X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
  828 FORMAT (21X,'S T R E S S E S   I N   I S O P A R A M E T R I C  ',&
             ' S O L I D   ( C I H E X',I2,2H ))
  829 FORMAT (2X,7HELEMENT,5X,4HGRID,11X,'STRESSES IN BASIC COORDINATE',&
             'SYSTEM',13X,12HDIR. COSINES)
  830 FORMAT (7X,2HID,4X,5HPOINT,8X,6HNORMAL,12X,5HSHEAR,10X,9HPRINCIPAL,&
             10X,1HA,4X,1HB,4X,1HC,4X,11HMEAN STRESS,5X,9HMAX SHEAR)
  831 FORMAT (13X,'C O M P L E X   S T R E S S E S   I N   I S O P A R',&
             ' A M E T R I C   S O L I D   ( C I H E X',I2,2H ))
  832 FORMAT (7X,2HID,3X,6HPOINTS,5X,8HNORMAL-X,9X,8HNORMAL-Y,9X,&
             8HNORMAL-Z,9X,8HSHEAR-XY,9X,8HSHEAR-YZ,9X,8HSHEAR-ZX)
  833 FORMAT (25X,'F O R C E   D I S T R I B U T I O N   I N   B A R',&
             '   E L E M E N T S,10X,11H( C B A R )')
  834 FORMAT (21H0    ELEMENT  STATION,9X,11HBEND-MOMENT,22X,&
             11HSHEAR FORCE,21X,5HAXIAL)
  835 FORMAT (7X,3HID.,5X,5H(PCT),5X,7HPLANE 1,8X,7HPLANE 2,11X,&
             7HPLANE 1,8X,7HPLANE 2,15X,5HFORCE,14X,6HTORQUE)
  836 FORMAT (25X,'S T R E S S   D I S T R I B U T I O N   I N   B A R',&
             '  E L E M E N T S,7X,11H( C B A R )')
  837 FORMAT (21H0    ELEMENT  STATION,4X,3HSXC,11X,3HSXD,11X,3HSXF,11X,&
              3HSXG,12X,5HAXIAL,10X,5HS-MAX, 9X,5HS-MIN,9X,4HM.S.)
  838 FORMAT (7X,3HID.,5X,5H(PCT))
  839 FORMAT (21X,'F O R C E S  F O R  T H E  Q U A D R I L A T E R A L', &
             '  T H I N  S H E L L     ( C Q U A D T S )')
  840 FORMAT (6X,2HEL,            36X,6HFORCES,51X,7HMOMENTS )
  841 FORMAT (6X,2HID,5X,5HPOINT,9X,2HFX,17X,2HFY,17X,2HFZ,17X,2HMX,17X,&
             2HMY,17X,2HMZ )
  842 FORMAT (17X,'F O R C E S   I N   T R I A N G U L A R   T H I N  ',&
             'S H E L L   E L E M E N T S   ( C T R S H L )')
  843 FORMAT (19X,'S T R E S S E S  F O R  T H E  Q U A D R I L A T E ',&
            'R A L  T H I N  S H E L L     ( C Q U A D T S )')
  844 FORMAT (3X,9HEL STRESS,8X,28HMEMBRANE  STRESS  RESULTANTS,24X,&
             17HFLEXURAL  MOMENTS,27X,5HSHEAR )
  845 FORMAT (3X,'ID  POINT   NORMAL(NX)     NORMAL(NY)     SHEAR(NXY)',&
             '     NORMAL(MX)     NORMAL(MY)     TORQUE(MXY)     ',&
             'NORMAL(QX)     NORMAL(QY)')
  846 FORMAT (18X,'S T R E S S E S   I N   T R I A N G U L A R   T H I',&
             ' N   S H E L L   E L E M E N T S   ( C T R S H L )')
  847 FORMAT (5X,'S T R E S S E S  I N  A X I S - S Y M M E T R I C  ',&
             'T R I A N G U L A R  R I N G  E L E M E N T S  (CTRIAAX)')
  848 FORMAT (' ELEMENT   HARMONIC    POINT    RADIAL      AXIAL',6X,&
             'CIRCUM.     SHEAR      SHEAR      SHEAR      F L U X   ',&
             'D E N S I T I E S', /,' ID.       NUMBER      ANGLE     ',&
             '(R)         (Z)     (THETA-T)    (ZR)       (RT)       ',&
             '(ZT)        (R)        (Z)        (T)')
  849 FORMAT (11X,'F O R C E S  I N  A X I S - S Y M M E T R I C  T R ',&
             'I A N G U L A R  R I N G  E L E M E N T S  (CTRIAAX)')
  850 FORMAT (1X,113H  ELEMENT   HARMONIC    POINT            RADIAL      &
              &      CIRCUMFERENTIAL            AXIAL                CHARGE,&
              /1X,'    ID.      NUMBER     ANGLE             (R)',17X,&
              '(THETA-T)                (Z)')
  851 FORMAT (5X,'S T R E S S E S  I N  A X I S - S Y M M E T R I C  T',&
            ' R A P E Z O I D A L  R I N G  E L E M E N T S  (CTRAPAX)')
  852 FORMAT (11X,'F O R C E S  I N  A X I S - S Y M M E T R I C  T R ',&
             'A P E Z O I D A L  R I N G  E L E M E N T S  (CTRAPAX)')
  853 FORMAT (43X,45HE L E M E N T   S T R A I N   E N E R G I E S )
  854 FORMAT (30X,15HELEMENT-TYPE = ,2A4,9X,23H* TOTAL FOR ALL TYPES = ,&
             1P,E16.7, /1H0,95X,1H*, /36X,10HELEMENT-ID,10X,&
             13HSTRAIN-ENERGY,11X,16HPERCENT OF TOTAL )
  855 FORMAT (42X,47HG R I D   P O I N T   F O R C E   B A L A N C E )
  856 FORMAT (11H   POINT-ID,4X,10HELEMENT-ID,5X,6HSOURCE,13X,2HT1,13X,&
             2HT2,13X,2HT3,13X,2HR1,13X,2HR2,13X,2HR3)
  857 FORMAT (22X,'F O R C E S   I N   T R I A N G U L A R   P L A T E',&
             '  E L E M E N T S   ( C T R P L T 1 )')
  858 FORMAT (20X,'S T R E S S E S   I N   T R I A N G U L A R   ',&
             'P L A T E   E L E M E N T S   ( C T R P L T 1 )')
  859 FORMAT (1H0,9X,7HELEMENT,4X,5HPOINT,7X,2(11HBEND-MOMENT, 9X),&
             12HTWIST-MOMENT,2(11X,5HSHEAR,4X))
  860 FORMAT (12X,3HID.,7X,3HNO.,13X,1HX,19X,1HY,39X,1HX,19X,1HY)
  861 FORMAT (1H0, 8H ELEMENT, 2X, 5HPOINT, 5X, 5HFIBER, 11X,&
             32HSTRESSES IN ELEMENT COORD SYSTEM, 12X,&
             31HPRINCIPAL STRESSES (ZERO SHEAR), 11X, 3HMAX)
  862 FORMAT (3X,3HID.,6X,3HNO.,5X,8HDISTANCE, 7X,8HNORMAL-X,6X,&
             8HNORMAL-Y,6X,8HSHEAR-XY,8X,5HANGLE,9X,5HMAJOR,9X,5HMINOR,&
             10X,5HSHEAR)
  863 FORMAT (18X,'S T R E S S E S   I N   T R I A N G U L A R   ',&
             'M E M B R A N E   E L E M E N T S   ( C T R I M 6 )')
  864 FORMAT (1H0, 8H ELEMENT, 5X, 5HPOINT, 7X,&
             32HSTRESSES IN ELEMENT COORD SYSTEM, 13X,&
             31HPRINCIPAL STRESSES (ZERO SHEAR), 13X, 3HMAX)
  865 FORMAT (4X,3HID.,8X,3HNO., 5X,8HNORMAL-X,7X,8HNORMAL-Y,7X,&
             8HSHEAR-XY,8X,5HANGLE,10X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
  866 FORMAT (2(24X, 6HMOMENT, 9X, 6HMOMENT, 15X), /2(6X, 7HSUBCASE,11X,&
             7HPTS 1,3, 8X, 7HPTS 2,4, 14X))
  867 FORMAT (6X, 7HSUBCASE, 2X, 2(11X, 11HBEND-MOMENT), 11X,&
             12HTWIST-MOMENT, 13X, 5HSHEAR, 17X, 5HSHEAR,&
             /31X, 1HX, 21X, 1HY, 43X, 1HX, 21X, 1HY)
  868 FORMAT (4(6X, 7HSUBCASE, 9X, 5HFORCE, 6X))
  869 FORMAT (5X, 7HSUBCASE, 11X, 3HSA1, 12X, 3HSA2, 12X, 3HSA3, 10X,&
             12HAXIAL-STRESS, 8X, 6HSA-MAX, 9X, 6HSA-MIN, 11X,6HM.S.-T,&
             /23X, 3HSB1, 12X, 3HSB2, 12X, 3HSB3, 30X, 6HSB-MAX, 9X,&
             6HSB-MIN, 11X, 6HM.S.-C)
  870 FORMAT (2(54X, 6HSAFETY), /2(5X, 7HSUBCASE, 14X, 7HMAXIMUM, 8X,&
             7HAVERAGE, 6X, 6HMARGIN))
  871 FORMAT (19X, 5HFIBRE, 11X, 32HSTRESSES IN ELEMENT COORD SYSTEM,&
             13X, 31HPRINCIPAL STRESSES (ZERO SHEAR), 10X, 7HMAXIMUM,/&
             5X, 7HSUBCASE, 6X, 8HDISTANCE, 7X, 8HNORMAL-X, 7X,&
             8HNORMAL-Y, 6X, 8HSHEAR-XY, 7X, 5HANGLE, 9X, 5HMAJOR,&
             11X, 5HMINOR, 10X, 5HSHEAR)
  872 FORMAT (4(6X, 7HSUBCASE, 8X, 6HSTRESS, 6X))
  873 FORMAT (107X, 22HOCTAHEDRAL    PRESSURE,/5X, 10HSUBCASE   , 8X,&
             8HSIGMA-XX, 6X, 8HSIGMA-YY, 6X, 8HSIGMA-ZZ, 7X, 6HTAU-YZ,&
             8X, 6HTAU-XZ, 8X, 6HTAU-XY, 8X, 5HTAU-0, 10X, 1HP)
  874 FORMAT (107X, 22HOCTAHEDRAL    PRESSURE,/5X, 11HSUBCASE    , 8X,&
             8HSIGMA-XX, 6X, 8HSIGMA-YY, 6X, 8HSIGMA-ZZ, 7X, 6HTAU-YZ,&
             8X, 6HTAU-XZ, 8X, 6HTAU-XY, 8X, 5HTAU-0, 10X, 1HP)
  875 FORMAT (32X,'F O R C E S   O F   M U L T I - P O I N T   C O N S',&
             ' T R A I N T')
  876 FORMAT (2X,7HELEMENT,4X,16HMAT. COORD. SYS.,6X,&
             33HSTRESSES IN MATERIAL COORD SYSTEM,12X,&
             31HPRINCIPAL STRESSES (ZERO SHEAR),12X,3HMAX)
  877 FORMAT (4X, 3HID., 6X, 15HID./OUTPUT CODE,&
             5X, 8HNORMAL-X, 7X, 8HNORMAL-Y, 6X, 8HSHEAR-XY,&
             7X, 5HANGLE, 9X, 5HMAJOR, 11X, 5HMINOR, 10X, 5HSHEAR)
  878 FORMAT (43X,'S T R E S S E S   A T   G R I D   P O I N T S')
  879 FORMAT (7X,'S T R A I N S / C U R V A T U R E S   I N   G E N E ',&
             'R A L   T R I A N G U L A R   E L E M E N T S',6X,&
             '( C T R I A 1 )')
  880 FORMAT (7X,'S T R A I N S / C U R V A T U R E S   I N   G E N E ',&
             'R A L   T R I A N G U L A R   E L E M E N T S',6X,&
             '( C T R I A 2 )')
 
END SUBROUTINE ofp1b
