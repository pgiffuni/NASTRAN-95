SUBROUTINE subbb
     
 COMPLEX :: sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,ai,a,b,bsycon,alp,f1,am1,aln,blkapm,bkdel3,  &
     f1s,c1,c2p,c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,  &
     sum2,f2,blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,  &
     gusamp,c1p,c1n,bkdel1,bkdel2,blkap1,arg,arg2,  &
     ft3tst,bc,bc2,bc3,bc4,bc5,ca1,ca2,ca3,ca4,clift,  &
     cmomt,pres1,pres2,pres3,pres4,qres4,cexp4c,fqa,  &
     fqb,t1,t2,t3,t4,cexp2a,cexp2b,cexp2c,cexp4a,  &
     cexp4b,fq7,c1a,c3a,c4a,const,cexp3,cexp4,cexp3a, cexp3b,cexp3c
 DIMENSION       pres1(21),pres2(21),pres3(21),pres4(21),qres4(21),  &
     sbkde1(201),sbkde2(201),sumsv1(201),sumsv2(201),  &
     svkl1(201),svkl2(201),xlsv1(21),xlsv2(21), xlsv3(21),xlsv4(21)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,ibbout
 COMMON /blk1  / scrk,sps,sns,dstr,ai,pi,del,sigma,beta,res
 COMMON /blk2  / bsycon
 COMMON /blk3  / sbkde1,sbkde2,f4,f4s,am4,f5s,f6s,am4tst,sum3,sum4,  &
     am5tt,am6,sumsv1,sumsv2,svkl1,svkl2,f5,f5t,am5,  &
     am5t,a,b,alp,f1,am1,aln,blkapm,bkdel3,f1s,c1,c2p,  &
     c2n,c2,amtest,ft2,blam1,ft3,am2,sum1,sum2,f2,  &
     blam2,ft2t,c1t,ft3t,f2p,am2p,sum1t,sum2t,c1p,c1n,  &
     bkdel1,bkdel2,blkap1,arg,arg2,ft3tst,bc,bc2,bc3,  &
     bc4,bc5,ca1,ca2,ca3,ca4,clift,cmomt,pres1,pres2,  &
     pres3,pres4,qres4,fqa,fqb,fq7
 COMMON /blk4  / i,r,y,a1,b1,c4,c5,gl,i6,i7,jl,nl,ri,rt,r5,sn,sp,  &
     xl,y1,amu,gam,idx,inx,nl2,rl1,rl2,rq1,rq2,xl1,  &
     alp1,alp2,gamn,gamp,iner,iout,redf,stag,step,  &
     amach,betnn,betnp,bkap1,xlsv1,xlsv2,xlsv3,xlsv4,  &
     alpamp,amoaxs,gusamp,disamp,pitaxs,pitcor
 
 s1    = 2.0 + sns - sps
 t1    = CEXP(-ai*sigma)
 t2    = CEXP(+ai*sigma)
 temp  = s1/rl2
 c1a   = ai*gl
 const = b*ai*(del-amu)*blam2/bkap1
 cexp3 = CEXP(c1a*sps )
 cexp4 = CEXP(c1a*temp)
 xl    = sps
 DO  jl = 1,nl2
   pres3(jl) = (ft2t*cexp3+ft3t+const*xl)*t1
   cexp3 = cexp3*cexp4
   xl    = xl + temp
 END DO
 ft3tst = 0.0
 ft2  = 0.0
 ft3  = 0.0
 ft2t = 0.0
 ft3t = 0.0
 fqa  = bkdel1/(bc*beta)*(a*ai*bkdel2/bkdel1-b*blkap1)*  &
     CEXP(-ai*(del*sps-sigma)/2.0)
 DO  i = 1,50
   rt  = 0.0
   r   = i - 1
   ri  = (-1.0)**(i-1)
!WKBR ALP = SQRT((R*PI/SNS)**2+SCRK**2)
   alp = SQRT((r*pi/sns)**2+scrk**2)
   aln = -alp
   CALL akapm (alp,bkdel3)
   t3  = alp - del
   svkl1(i) = bkdel3
   IF (i == 1) rt = 1.0
   sum1  = (alp-amu)/(t3)*(ri-CEXP(ai*(t3)*sps)*t2)/  &
       (beta*(1.0+rt))*ri/(sns*alp)*bkdel1/bkdel3*(a*ai*bkdel2/  &
       bkdel1*(t3)/(t3+gl)-b*blkap1-b/(t3))
   sum1t = (alp-amu)/(t3)*(1.0-CEXP(ai*(t3)*sps)*t2*ri)/  &
       (beta*(1.0+rt))*ri/(sns*alp)*bkdel1/bkdel3*(a*ai*bkdel2/  &
       bkdel1*(t3)/(t3+gl)-b*blkap1-b/(t3))
   sumsv1(i) = (alp-amu)/(t3)*(1.0-CCOS((t3)*sps+sigma+r*pi))/  &
       (beta*(1.0+rt)*sns*alp)*bkdel1/bkdel3*CEXP(-2.0*ai*(alp-  &
       del))*(a*bkdel2/bkdel1*(t3)/(t3+gl)+b*ai*blkap1+b*ai/(t3))
   ft2   = sum1*ai/(t3)*(CEXP(-2.0*ai*(t3))-CEXP(-ai*(sps-sns)*(t3))) + ft2
   ft3   = sum1*(2.0*ai*CEXP(-2.0*ai*(t3))/(t3)-ai*(sps-sns)/  &
       (t3)*CEXP(-ai*(t3)*(sps-sns))+CEXP(-2.0*ai*(t3))/  &
       ((t3)**2)-CEXP(-ai*(t3)*(sps-sns))/((t3)**2)) + ft3
   ft2t  = sum1t*t1*CEXP(-ai*(t3)*sps)*ai/(t3)*(CEXP(-ai*(t3)*(s1))-  &
       1.0) + ft2t
   ft3t  = sum1t*t1*CEXP(-ai*(t3)*sps)*((s1)*ai/(t3)*CEXP(-ai*(t3)*  &
       (s1)) + 1.0/((t3)**2)*(CEXP(-ai*(t3)*(s1))-1.0)) + ft3t
   CALL akapm (aln,bkdel3)
   t4    = aln - del
   svkl2(i) = bkdel3
   sum2  = (aln-amu)/(t4)*(ri-CEXP(ai*(t4)*sps)*t2)/(beta*(1.0+rt))*  &
       ri/(sns*aln)*bkdel1/bkdel3*(a*ai*bkdel2/bkdel1*(t4)/  &
       (t4+gl)-b*blkap1-b/(t4))
   sum2t = (aln-amu)/(t4)*(1.0-CEXP(ai*(t4)*sps)*t2*ri)/(beta*(1.0+  &
       rt))*ri/(sns*aln)*bkdel1/bkdel3*(a*ai*bkdel2/bkdel1*(t4)/  &
       (t4+gl)-b*blkap1-b/(t4))
   sumsv2(i) = (aln-amu)/(t4)*(1.0-CCOS((t4)*sps+sigma+r*pi))/  &
       (beta*(1.0+rt)*sns*aln)*bkdel1/bkdel3*CEXP(-2.0*ai*(t4))*  &
       (a*bkdel2/bkdel1*(t4)/(t4+gl)+b*ai*blkap1+b*ai/(t4))
   ft2   = ft2+sum2*ai/(t4)*(CEXP(-2.0*ai*(t4))-CEXP(-ai*(sps-sns)* (t4)))
   ft2t  = sum2t*t1*CEXP(-ai*(t4)*sps)*ai/(t4)*(CEXP(-ai*(t4)*(s1))-  &
       1.0) + ft2t
   ft3   = ft3+sum2*(2.0*ai*CEXP(-2.0*ai*(t4))/(t4)-ai*(sps-sns)/  &
       (t4)*CEXP(-ai*(t4)*(sps-sns))+CEXP(-2.0*ai*(t4))/  &
       ((t4)**2)-CEXP(-ai*(t4)*(sps-sns))/((t4)**2))
   ft3t  = ft3t+sum2t*t1*CEXP(-ai*(t4)*sps)*((s1)*ai/(t4)*  &
       CEXP(-ai*(t4)*(s1))+1./((t4)**2)*(CEXP(-ai*(t4)*(s1))-1.))
   i7    = i
   aa    = sps - sns
   temp  = s1/rl2
   temp2 = r*pi/sns
   const = 4.0/pi*fqa
   temp3 = r + rt
   c3a   = -ai*t3
   c4a   = -ai*t4
   c1a   = ai*del
   cexp3a = CEXP(c3a*aa)
   cexp3b = CEXP(c3a*sps)
   cexp3c = CEXP(c3a*temp)
   cexp4a = CEXP(c4a*aa)
   cexp4b = CEXP(c4a*sps)
   cexp4c = CEXP(c4a*temp)
   cexp2a = CEXP(c1a*aa)
   cexp2b = CEXP(c1a*sps)
   cexp2c = CEXP(c1a*temp)
   xl1    = aa
   DO  jl = 1,nl2
     pres2(jl) = sum1*cexp3a+sum2*cexp4a + pres2(jl)
     pres2(jl) = pres2(jl) + const*cexp2a*ri/temp3*SIN(temp2*(xl1-sps))
     xl2 = xl1 + sns
     pres3(jl) = (sum1t*cexp3b+sum2t*cexp4b)*t1 + pres3(jl)
     pres3(jl) = pres3(jl)+const*cexp2b/temp3*SIN(temp2*(xl2-sps))*t1
     xl1 = xl1 + temp
     cexp3a = cexp3a*cexp3c
     cexp4a = cexp4a*cexp4c
     cexp2a = cexp2a*cexp2c
     cexp3b = cexp3b*cexp3c
     cexp4b = cexp4b*cexp4c
     cexp2b = cexp2b*cexp2c
   END DO
   IF (cabs((ft3-ft3tst)/ft3) < 0.0006) GO TO 65
   ft3tst = ft3
 END DO
 GO TO 9994
 65 CONTINUE
 ft3tst = ft3
 f2     = f2  + ft2
 am2    = am2 + ft3
 f2p    = f2p + ft2t
 am2p   = am2p+ ft3t
 aa     = sps - sns
 aa1    = sps + sns
 aa2    = sps + 2.0*sns
 temp   = s1/rl2
 xl     = aa
 c1a    = ai*del
 cexp3  = CEXP(c1a*aa)
 cexp3c = CEXP(c1a*temp)
 cexp4  = CEXP(c1a*sps)
 const  = 2.0*fqa
 cexp2a = t1*const
 DO  jl = 1,nl2
   step = 0.0
   IF (xl >= aa1) step = 1.0
   pres2(jl) = pres2(jl) + const*cexp3*((xl-sps)/sns-2.0*step)
   xl2  = xl + sns
   step = 0.0
   IF (xl2 >= aa2) step = 1.0
   pres3(jl) = pres3(jl) - cexp2a*cexp4*(1.0-(xl2-sps)/sns+2.0*step)
   cexp3 = cexp3*cexp3c
   cexp4 = cexp4*cexp3c
   xl  = xl + temp
 END DO
 gam = sps*del - sigma
 c1p = (gam/dstr) - scrk
 c2p = (gam/dstr) + scrk
 alp = gam*sps/(dstr**2) - sns/dstr*CSQRT(c1p)*CSQRT(c2p)
 t3  = alp - del
 f4  = CEXP(ai*(alp*sps-gam))*(alp*sps-gam)/((alp*dstr**2-gam*sps)* (t3))
 CALL akapm (alp,bkdel3)
 sbkde1(1) = bkdel3
 sbkde2(1) = 0.0
 CALL akappa (del,bkap1)
 carg = del - gl
 CALL akappa (carg,ckap1)
 f4  = f4*bkdel3/(bkdel1*bkap1)*(a*(bkdel1/bkdel2*(t3)/(t3+gl)*  &
     (del-gl-amu)*CEXP(2.0*ai*gl)*bkap1/ckap1)+b*ai*(1.0-2.0*ai*  &
     (del-amu)-(del-amu)*res)-b*ai*(del-amu)*(blkap1-1.0/(t3)))
 f5s = b*ai/(bkdel1*bkap1)*(1.0-2.0*ai*(del-amu) - (del-amu)*res -  &
     (del-amu)*blkap1)
 f6s = a/(bkdel1*bkap1)*(bkdel1/bkdel2*(del-gl-amu)*CEXP(2.0*ai*gl)  &
     *bkap1/ckap1)
 f4s = f4
 fq7 = bc*(f6s+f5s)
 temp  = (sps-sns)/rl1
 temp2 = 2.0 - sps
 const = -t1*f4s
 c1a   = -ai*t3
 cexp3a = CEXP(c1a*sns)
 cexp3b = CEXP(c1a*temp)
 DO  jl = 1,nl
   pres4(jl) = const*cexp3a
   cexp3a = cexp3a*cexp3b
 END DO
 c1  = CEXP(-ai*(t3)*sps)
 c2  = CEXP(-ai*(t3)*sns)
 f4  = f4*ai*t1/(t3)*(c1-c2)
 am4 = f4s*t1*(ai*sps*c1/(t3)-ai*sns*c2/(t3)+(c1-c2)/  &
     ((t3)**2))+f4s*ai*(2.0-sps)*t1/(t3)*(c1-c2)
 CALL subc
 RETURN
 
 9994 WRITE  (ibbout,3015) ufm
 3015 FORMAT (a23,' - AMG MODULE -SUBROUTINE SUBC.  AM4 LOOP DID NOT ',  &
     'CONVERGE.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE subbb
