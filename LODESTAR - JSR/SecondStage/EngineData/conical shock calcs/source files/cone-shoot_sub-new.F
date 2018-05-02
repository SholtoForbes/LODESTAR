      subroutine cone_shoot(M1,thetac,alpha,Mc,pcop1,TcoT1)
c
c    M1 - flight Mach number
c    thetac - cone half angle (degrees)
c    alpha - angle-of-attack (degrees)
c
c     PROGRAM cone-shoot
c
c   program to calculate the properties on a cone
c   at a specified supersonic Mach number
c
c   Uses a Newton Raphson solver to choose a
c   starting shock angle for the integration of the
c   Taylor-Macoll equations on a cone of given
c   half angle.
c   (1) cone-shoot does the shooting for the initial shock angle
c   (2) funcv performs the integration from a given shock angle to the cone half angle.
c   (3) answer is when the angular velocity at the surface the cone is smaller than some tolerance.
c
c   Author: M.Smart
c   Date: 8 October 2016
c
c   Notes:  18 October 2016 - try using Brents method instead of Newt, as know
c                             that cone angle will be in between 2-D shock angle and
c                             Mach angle.
c                           - This works better than original method using newt.
c                           - need factor of 1.01 on Mach angle, otherwise sometimes blows up
c
c  uses:  oshock, brent, newt
c
      INTEGER i,m,n
      REAL dx,x1,x2
      real M1,Mach,Mc
c
      external cone
c
      common /requirement/ Mach,hangle,vrc,pt2opt1
      common /gas/ gam,g1,g2
c
      pi = acos(-1.0000000)
c
      gam = 1.400000000
      g1 = gam + 1.0000000
      g2 = gam - 1.0000000
c
      Mach = M1
c
      dx=1.e-4
c
      hangle = thetac + alpha
c
      write(6,'(a,f8.2,a,f5.2)')'  M1=',Mach,
     $             '  Equiv. cone half angle=',hangle
c
c  2-D Wedge Shock angle always greater than conical shock angle
c  Use this as upper limit for shock angle.  Use Mach angle as lower
c  limit for cone shock angle
c
      theta = hangle
      call oshock(M1,theta,gam,beta,mprat,vrat,trat,M2)
      thetasmax = beta*pi/180.00000
      thetasmin = 1.01*asin(1.000000/M1)
c
      tol = 1.0e-8
      thetas = brent(cone,thetasmin,thetasmax,tol)
c                              (deg.)
        Mc = (g2/2.000000*((1.000000/Vrc)**2 - 1.000000))**(-0.500000)
        pcopt2 = (1.00000 + 0.500000*g2*Mc**2)**(-gam/g2)
        pt1op1 = (1.00000 + 0.500000*g2*M1**2)**(gam/g2)
        pcop1 = pcopt2*pt2opt1*pt1op1
        TcoTt1 = 1.000000/(1.00000 + 0.500000*g2*Mc**2)
        Tt1oT1 = (1.00000 + 0.500000*g2*M1**2)
        TcoT1 = TcoTt1*Tt1oT1
        write(6,*) 'Cone: Vrat,Mach,pc/p1,Tc/T1',Vrc,Mc,pcop1,TcoT1
c
      END
c
c
c
      subroutine oshock(M1,theta,gam,beta,mprat,vrat,trat,M2)
c
*********************************
c  Oblique SHOCK.f **************
c *******************************
c
c  program to calculate the weak oblique shock properties given
c  M, theta (turning angle, deg.) and gamma.
c
c  author - M.Smart
c  date - 18 January 2001
c
c  Notes - 22 September 2106 - add beta as an output
c
      real M,M1,M2
c
      external shock
c
      common /os/ M,th,gamm
c
      g1 = gam + 1.000000000
      g2 = gam - 1.000000000
      pi = acos(-1.000000000)
      M = M1
      th = theta
      gamm = gam
c
c      write(6,*) M1,theta,gam
      bmin = 0.10000000 +
     $       asin(1.00000000/M1)*180.000000/pi
      bmax = 60.000000
      tol = 1.0e-8
      beta = brent(shock,bmin,bmax,tol)
c                              (deg.)
c      write(6,*)' beta = ',beta
      t1 = sin((beta-theta)*pi/180.0000)**2
      t2 = g2*(M1*sin(beta*pi/180.000000))**2 + 2.00000000
      t3 = 2.00000000*gam*(M1*sin(beta*pi/180.000000))**2 - g2
      t4 = g1*(M1*sin(beta*pi/180.000000))**2
      t5 = (M1*sin(beta*pi/180.000000))**2 - 1.0000000
      t6 = gam*(M1*sin(beta*pi/180.000000))**2 + 1.00000000
      M2 = sqrt(t2/t3/t1)
      prat = t3/g1
      ptrat = (t4/t2)**(gam/g2)*(g1/t3)**(1.0000000/g2)
      trat = t3*t2/g1/t4
      vrat = sqrt(1.000000-4.000000*t5*t6/(t4*g1*M1**2))
c      write(6,*) '     M1      theta     gamma'
c      write(6,'(3f10.4)')M1,theta,gam
c      write(6,*) '   beta       M2       p2/p1   pt2/pt1   t2/t1'
c      write(6,'(4f10.4)')beta,M2,prat,ptrat,trat
c
      return
      end
c
c
c
      FUNCTION brent(func,x1,x2,tol)
c
      INTEGER ITMAX
c
      real brent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      real a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
        write(6,*)'roots not bracketed'
        write(6,*)'fa,fb:',fa,fb
c        pause
c     *'root not bracketed by initial limits'
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.50000000*tol
        xm=.5000000*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.0000000)then
          brent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.0000000*xm*s
            q=1.0000000-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.00000000*xm*q*(q-r)-(b-a)*(r-1.0000000))
            q=(q-1.000000)*(r-1.000000)*(s-1.0000000)
          endif
          if(p.gt.0.0000000000) q=-q
          p=abs(p)
          if(2.00000*p .lt. min(3.0000*xm*q-(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'brent exceeding maximum iterations'
      brent=b
      return
      END
c
c
c
      function shock(x)
c
      real M
c
      common /os/ M,th,gamm
c
      pi = acos(-1.0000000000)
      g1 = gamm + 1.000000000
      g2 = gamm - 1.000000000
      beta = x
c
      t1 = g1*M**2
      t2 = 2.0000000*
     $     ((M*sin(beta*pi/180.0000000))**2-1.000000)
      t3 = 1.0000000/tan(th*pi/180.000000)
      shock = t3 - tan(beta*pi/180.000000)*(t1/t2 - 1.0000000)
c      write(6,*)'M1,theta,beta,shock:',M,th,beta,shock
c
      return
      end
c
c
c
      function cone(thetasr)
c
      INTEGER n2,nvar,kmax,kount,KMAXX,NMAX
      REAL x1,x2,dxsav,EPS
      PARAMETER (n2=1,nvar=2,NMAX=50,KMAXX=2000,EPS=1.e-6)
      REAL h1,hmin,ystart(nvar),xp(KMAXX),yp(NMAX,KMAXX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
c
CU    USES derivs,odeint,rkqs
c
      INTEGER nbad,nok
      EXTERNAL derivs,rkqs
c
      real Mach,M2
c
      common /requirement/ Mach,hangle,vrc,pt2opt1
      common /gas/ gam,g1,g2
c
      pi = acos(-1.0000000000)
c
c  perform integration; note that starting point of the integration
c  is the parameter we are trying to find.
c
      x1=thetasr
      x2=hangle*pi/180.00000
      thetas = thetasr*180.0000/pi
c
      call oshock_prop(Mach,thetas,M2,delta,vr2,vtheta2,pt2opt1)
      ystart(1)=vr2
      ystart(2)=-vtheta2
c      write(6,*) Mach,thetas,delta,M2,vr2,vtheta2
      h1=(x2-x1)/100.
      hmin=0.000000
      kmax=KMAXX
      dxsav=min((x2-x1)/100.,0.1)
      call odeint(ystart,NVAR,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
c      write(6,'(/1x,a,t30,i8)') 'Successful steps:',nok
c      write(6,'(1x,a,t30,i8)') 'Bad steps:',nbad
c      write(6,'(1x,a,t30,i3)') 'Function evaluations:',nrhs
c      write(6,'(1x,a,t30,i8)') 'Stored intermediate values:',kount
c
c      do 100 i = 1,kount
c        write(6,*) xp(i)*180.0000/pi,yp(1,i),yp(2,i)
c  100 continue
      vrc = yp(1,kount)
      vthetac = -yp(2,kount)
c                          vthetac = -yp2 ! (Justine)
c
c  vtheta2 = zero at cone surface, so descrepancy, f, is vtheta2
c
      write(6,*)' Shock angle, Vtheta on cone:',x1*180./pi,vthetac
      cone = vthetac
c
      return
c
      END
c
c
c
      subroutine oshock_prop(M1,thetas,M2,delta,vr2,vtheta2,pt2opt1)
c
c   Calculate the flow properties behind a shock using oblique shock relations
c   given M1 and the shock angle
c
      real M1,Mn1,Mn2,M2
c
      common /gas/ gam,g1,g2
c
      pi = acos(-1.00000000)
c
      thetasr = thetas*pi/180.0000000
c
      Mn1=M1*sin(thetasr)
      t2=2.00000000+g2*(Mn1**2)
      t3=2.000000*gam*(Mn1**2)-g2
      t4=g1*(Mn1**2)
      Mn2=(t2/t3)**(0.50000)
      deltar=atan((2.000000/(tan(thetasr)))*((Mn1**2.0000)-1.000000)/
     $         (((M1**2)*gam+(M1**2)*cos(2.000000*thetasr)+2.000000)))
      M2=Mn2/(sin(thetasr-deltar))
      V2=((2.000000/(g2*(M2**2)))+1.00000)**(-0.500000)
      vr2=V2*cos(thetasr-deltar)
      vtheta2=V2*sin(thetasr-deltar)
      pt2opt1=(g1/t3)**(1.0000/g2)*
     $        (g1*Mn1**2/t2)**(gam/g2)
c
      delta = deltar*180.00000/pi
c
      return
c
      end
c
c
c
      SUBROUTINE derivs(x,y,dydx)
c
      REAL x,y(*),dydx(*)
c
      common /gas/ gam,g1,g2
c
      pi = acos(-1.0000000)
c
c  ode's:
c            x - theta (radians)
c            y1 = normalised radial velocity component
c            y2 = dy1/dx = normalised angular velocity component (irrotational condition)
c
      dydx(1)= y(2)
      t1 = 1.0000000 - y(1)**2 - y(2)**2
      cot = 1.000000/(tan(x))
      t2 = 2.000000*y(1) + y(2)*cot
      t3 = -g2/2.0000000*t1*t2 + y(2)**2*y(1)
      t4 = g2/2.0000000*t1 - y(2)**2
c      t3 = -g2*t1*t2 + y(2)**2*y(1)
c      t4 = g2*t1 - y(2)**2
      dydx(2)= t3/t4
c      write(6,*) x*180.00/pi,y(1),y(2),dydx(2)
c
      return
      END
c
c
c
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=2000,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
c
c   pauses commented out
c
c        if(abs(hnext).lt.hmin) pause
c     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
c      pause 'too many steps in odeint'
      return
      END
C
C
C
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      REAL errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,
     *ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x)then
c        pause 'stepsize underflow in rkqs'
        goto 1
        endif
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
C
C
C
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
