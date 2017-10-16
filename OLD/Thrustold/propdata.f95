cc
      subroutine thrust_5_int(M1,p1,T1,wcap,fs,f,isp,phi,wf,
     $                        Mache,prese,tempe,gammae,Re)
c
c  Input - M1 - Mach number entering engine
c          p1 - pressure intering engine (Pa)
c          T1 - Temperature entering engine (K)
c          wcap - capture width of engine (m)
c
c   Output - fs - specific thrust (m/s)
c            f - thrust (N)
c            isp - specific impusle (s)
c            phi - equivalence ratio
c            wf - fuel mass flow (kg/s)
c            Mache - nozzle exit Mach number
c            prese - nozzle exit pressure (Pa)
c            Tempe - nozzle exit Temperature (K)
c            gammae - nozzle exit ratio of specific heats 
c            Re - nozzle exit gas constant (J/kgK)
c
c
c   Free format code to be compatable with f90 interpolation routine
c   Lahey requires that all code linked together be either fixed-format (.f or .for postscript)
c   or free format (.f90 or .f95 postrcript).
c   Simply just replace c with ! in first column of each comment line
c
c   Date 1 March 2010
c   Author - M. Smart
c
c  notes:  25 October - introduce more cycle data to improve interpolation (particularly phi)
c
c   OLD NOTES
c *********************************************************************************
c
c  date - 23 February 2010 (thrust5_int)
c
c   - remove dependence of engine performance on M and alpha.  Replace with 
c     properties entering the engine - M1 p1, and T1. M1 and T1 togther define 
c     the total and local flow conditions entering the engine; these dictate 
c     the performance of the engine.  However, the engine performance is 
c     calculated assuming a value of p1 (based on q = 50 kpa, thetab = 6 deg.
c     and an angle of attack), so actual p1 is needed to scale the exit pressure
c     of the engine (Prese).
c
c   - all input based on flight conditions has been removed (change to input file;
c     now thrust5.in
c   - forebody calcs removed.
c***********************************************************************************
c  
c
      real M1,M2,mc,isp,M1des,Mrati,Mache
C Dawid Plotting 15/10/2011
      common /ispplot/ w2
c
c
c
c
      open(888, file = 'dawidmpt1')
      write(888,*)m1, p1, t1
c
!c     RESTM12 engine
!c
!      wcapstandard = 0.200
!c                     meters
!c
!      R0 = 0.1774982*wcap/wcapstandard
!c                     meters
!c
!      Acap = 0.635005*R0**2
!c                      m^2

!      FULL CAPTURE ENGINE
      wcapstandard = 0.2156
!                       meters
      R0 = wcap**2/wcapstandard**2
c                     meters
c
      Acap = 0.0470*R0
c                      m^2          


      M1des = 9.111
c
      if(M1 .gt. M1des) then
        write(6,*)' Inlet Over-sped'
        call inlet_restm12(M1des,ptrati,mc,Mrati,prati,Trati)
      else
        call inlet_restm12(M1,ptrati,mc,Mrati,prati,Trati)
      endif
c
      gam0 = 1.4000000
      r = 287.035
c                   J/kgK
c
c      w1 = mc*Acap*p1*M1*sqrt(gam0/r/T1)*3.0*1.26 
cc                      3 engine modules, full capture +26% DP10/2/2012
       w1 = mc*Acap*p1*M1*sqrt(gam0/r/T1)*4.0 
c        4 engine modules in Full capture engine 
C                      (kg/s)
c
      p2 = p1*prati
      T2 = T1*trati
      M2 = M1*Mrati
      w2 = w1
c
CCC      WRITE(*,*)'sendtorest',M1, T1
      call restm12(M1,T1,p1nom,isp,fs,phi,Mache,prese,tempe,gammae,Re)
c
c   (M1, and T1 are input variables  - 23 February 2010)
c
c    Must scale prese by ratio by P1/P1nom (assume P1 has no other effect on engine performance)
c
      prese = prese*p1/p1nom
c
      f = fs*w2
c            (N)
c
c  H2 fuel 
c
      fst = 0.0291
      wf = fst*w2*phi
c             (kg/s)
c
CCC      WRITE(*,*)'thrust',M1,f
c     dawid exit conditions 4/2012 
c      write (*,*) Mache,prese,tempe,gammae
c
c     NOZZLE
c      write(*,*)'thrust', mache

      return
      end
c
c
c
      subroutine restm12(xi,yi,p1nom,isp,fs,phi,
     $                   Mache,prese,tempe,gammae,Re)
c
c  fortran 90 subroutine to calculate the internal performance of the RESTM12 
c  flowpath (given M1 and T1) based on curve fits through a database.
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c **********************************************************************
c
c  xi: M1
c  yi: T1 (K)
c
c data - 1 March 2010
c author - M. Smart
c
c  notes:  20 October 2011 - add M7 al = -4,-2,0,2,4
c          24 October 2011 - simplfiy phi interpolation by only using
c                            interpolation below the Mach 7 line (second
c                            bivar subroutine: bivar_phi)
c       :  31 October 2011 - separate into two bivars; one below M = 7 and one above
c
      integer, parameter :: ni = 1
      real xi(ni),yi(ni),e1(ni),e2(ni),e3(ni),e4(ni),e5(ni)
      real e6(ni),e7(ni),e8(ni),e9(ni)
      real isp,Mache
C Dawid print time, 7/11/2011
C      COMMON C(3510)      
C      EQUIVALENCE (C(2000), TIME)
C      EQUIVALENCE (C(1203),ALPHAX)
c
      AL =  2.2235906196700304E+002
      BL =  2.4099186534551662E+000
      CL =  9.7131398721484852E-002
      DL = -7.4022839183780006E-004
      tlow  = AL + BL*xi(1) + CL*xi(1)**2 + DL*xi(1)**3
c
      Ah = -5.0131819173867848E+002
      Bh =  4.2638833081252193E+002
      Ch = -7.7515246232795320E+001
      Dh =  5.3058381342009682E+000
      thigh  = ah + bh*xi(1) + ch*xi(1)**2 + dh*xi(1)**3	  
c
      A6 =  1.8971714805999188E+003
      B6 = -6.3913717005006777E+002
      C6 =  8.4144214778663596E+001
      D6 = -4.0355209165118149E+000
      tlowM6 = a6 + b6*xi(1) + c6*xi(1)**2 + d6*xi(1)**3  
c
      A7 =  1.5056087084349767E+003
      B7 = -3.1406783229349946E+002
      C7 =  1.7498745858638578E+001
      D7 =  1.7781746014793676E-001
      tlowM7 = a7 + b7*xi(1) + c7*xi(1)**2 + d7*xi(1)**3
c
      A8 =  2.3916744281917681E+003
      B8 = -6.6295434011836960E+002
      C8 =  7.0588698211823271E+001
      D8 = -2.7076685693339320E+000
      tlowM8 = a8 + b8*xi(1) + c8*xi(1)**2 + d8*xi(1)**3
c
      A10 =  2.8743890187560614E+003
      B10 = -6.8809190345325840E+002
      C10 =  6.2664485135644767E+001
      D10 = -2.0420230248445130E+000
      tlowM10 = a10 + b10*xi(1) + c10*xi(1)**2 + d10*xi(1)**3
c
      A12 =  3.3145474496092438E+003
      B12 = -7.0132929401983256E+002
      C12 =  5.6076770367764112E+001
      D12 = -1.5966010571734792E+000
      tlowM12 = a12 + b12*xi(1) + c12*xi(1)**2 + d12*xi(1)**3
c
      A14 =  3.8296210675172551E+003
      B14 = -7.3814732153015211E+002
      C14 =  5.3380071137661240E+001
      D14 = -1.3664310805553339E+000
      tlowM14 = a14 + b14*xi(1) + c14*xi(1)**2 + d14*xi(1)**3
c	  
c      if((yi(1) .gt. thigh).or.(yi(1) .lt. tlow)) then
c        write(6,*)'M1=',xi(1),' out of range of propulsion deck'
C        write(6,*)'time, AoA', time, alphax
C        write(6,*)'thigh, tlow, yi(1)', thigh, tlow, yi(1)		
C        stop
c      elseif(yi(1) .lt. tlowM6) then
      if(yi(1) .lt. tlowM6) then
        write(6,*)'M1=',xi(1),' below range of propulsion deck'		
        ! IF THIS COMES UP DURING A TRAJECTORY with M1=0:
        !   YOU forgot to compile in dbl precision, use -dbl in LF95!!
        stop
        call bivar_67(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)		
      elseif((yi(1) .lt. tlowM7) .and. (xi(1) .lt. 7.0)) then
C          write(6,*)'M1=',xi(1),'less than M7'
C          write(6,*)'time', time
          call bivar_67(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)	
      elseif(yi(1) .lt. tlowM8) then
C          write(6,*)'M1=',xi(1),'less than M8'	  
          call bivar_78(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tlowM10) then
C          write(6,*)'M1=',xi(1),'less than M10'	  
          call bivar_810(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tlowM12) then
C          write(6,*)'M1=',xi(1),'less than M12'	  
          call bivar_1012(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tlowM14) then
C          write(6,*)'M1=',xi(1),'less than M14'	  
          call bivar_1214(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)	
      else
        write(6,*)'M1=',xi(1),' above range of propulsion deck'		
        stop
      endif		
c
      p1nom = e1(1)
      fs = e2(1)
      isp = e3(1)
      phi = e4(1)
      Mache = e5(1)
      prese = e6(1)
      tempe = e7(1)
      gammae = e8(1)
      Re = e9(1)
c
c     NOZZLE
c      write(*,*)'inletM12', mache

      return
      end
c
c
c
      subroutine 
     $   bivar_67(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c date - 5 November 2011
c author - M. Smart
c
      integer, parameter :: ni = 1
      integer, parameter :: nd = 10
      integer, parameter :: niwk = 31 * nd + ni
      integer, parameter :: nwk = 8 * nd
c
      integer i
      integer iwk(niwk)
      integer md
      real wk(nwk)
      real xd(nd)
      real xi(ni)
      real yd(nd)
      real yi(ni)
c
      real fsa(nd),ispa(nd),er(nd),Ma(nd),pa(nd),Ta(nd)
      real ga(nd),ra(nd),p1noma(nd)
      real fs(ni),isp(ni),phi(ni),Mache(ni),prese(ni)
      real Tempe(ni),gammae(ni),re(ni),p1nom(ni)
c
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M6; (6-10)M7
c
      xd(1) = 5.719		
      xd(2) = 5.449		
      xd(3) = 5.182	
      xd(4) = 4.915		
      xd(5) = 4.648		
      xd(6) = 6.663		
      xd(7) = 6.281		
      xd(8) = 5.931		
      xd(9) = 5.582	
      xd(10) = 5.234		  
c
      yd(1) = 239.20	
      yd(2) = 259.98	
      yd(3) = 283.15	
      yd(4) = 309.35	
      yd(5) = 339.08	
      yd(6) = 242.46	
      yd(7) = 267.27	
      yd(8) = 295.65	
      yd(9) = 328.56	
      yd(10) = 366.67	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 2.653e+3   
      p1noma(2) = 3.495e+3  
      p1noma(3) = 4.535e+3  
      p1noma(4) = 5.791e+3  
      p1noma(5) = 7.278e+3  
      p1noma(6) = 2.040e+3  
      p1noma(7) = 2.801e+3    
      p1noma(8) = 3.766e+3    
      p1noma(9) = 4.956e+3    
      p1noma(10) = 6.386e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 664.1
      fsa(2) = 617.9
      fsa(3) = 569.2
      fsa(4) = 541.7
      fsa(5) = 490.9
      fsa(6) = 550.6
      fsa(7) = 548.1
      fsa(8) = 557.6
      fsa(9) = 569.1
      fsa(10) = 583.3	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 2322.0
      ispa(2) = 2323.0
      ispa(3) = 2341.9
      ispa(4) = 2367.2
      ispa(5) = 2382.2
      ispa(6) = 1925.4	  
      ispa(7) = 1917.1
      ispa(8) = 1949.8
      ispa(9) = 1990.0
      ispa(10) = 2040.0	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      er(1) = 1.000
      er(2) = 0.93
      er(3) = 0.85
      er(4) = 0.80
      er(5) = 0.721
      er(6) = 1.0
      er(7) = 1.0
      er(8) = 1.0
      er(9) = 1.0
      er(10) = 1.0 
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, er, ni, xi, yi, phi, iwk, wk )
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 2.804
      Ma(2) = 2.809
      Ma(3) = 2.823
      Ma(4) = 2.830
      Ma(5) = 2.839
      Ma(6) = 2.959
      Ma(7) = 2.914
      Ma(8) = 2.916
      Ma(9) = 2.918
      Ma(10) = 2.919	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
      
c     NOZZLE
c      write(*,*)'bivar67', mache
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 10677.
      pa(2) = 12227.
      pa(3) = 13612.
      pa(4) = 15049.
      pa(5) = 14818.
      pa(6) = 9802.
      pa(7) = 12141.
      pa(8) = 14254.
      pa(9) = 16270.
      pa(10) = 18067.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1348.0
      ta(2) = 1302.7
      ta(3) = 1250.1
      ta(4) = 1217.2
      ta(5) = 1085.8
      ta(6) = 1427.9
      ta(7) = 1454.8
      ta(8) = 1452.5
      ta(9) = 1451.3
      ta(10) = 1451.7	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.2880
      ga(2) = 1.2922
      ga(3) = 1.2970
      ga(4) = 1.300
      ga(5) = 1.312
      ga(6) = 1.2826	
      ga(7) = 1.281
      ga(8) = 1.281
      ga(9) = 1.281
      ga(10) = 1.281	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit 
c
      ra(1) = 351.61
      ra(2) = 347.50
      ra(3) = 342.53
      ra(4) = 339.51
      ra(5) = 327.08
      ra(6) = 350.44
      ra(7) = 350.42
      ra(8) = 350.43
      ra(9) = 350.44
      ra(10) = 350.42	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ra, ni, xi, yi, re, iwk, wk )
c
      return
      end	
c
c
c
      subroutine 
     $  bivar_78(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c date - 5 November 2011
c author - M. Smart
c
      integer, parameter :: ni = 1
      integer, parameter :: nd = 10
      integer, parameter :: niwk = 31 * nd + ni
      integer, parameter :: nwk = 8 * nd
c
      integer i
      integer iwk(niwk)
      integer md
      real wk(nwk)
      real xd(nd)
      real xi(ni)
      real yd(nd)
      real yi(ni)
c
      real fsa(nd),ispa(nd),er(nd),Ma(nd),pa(nd),Ta(nd)
      real ga(nd),ra(nd),p1noma(nd)
      real fs(ni),isp(ni),phi(ni),Mache(ni),prese(ni)
      real Tempe(ni),gammae(ni),re(ni),p1nom(ni)
c
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M7; (6-10)M8
c
      xd(1) = 6.663		
      xd(2) = 6.281		
      xd(3) = 5.931	
      xd(4) = 5.582		
      xd(5) = 5.234		
      xd(6) = 7.534		
      xd(7) = 7.078		
      xd(8) = 6.642		
      xd(9) = 6.199	
      xd(10) = 5.762		  
c
      yd(1) = 242.46	
      yd(2) = 267.27	
      yd(3) = 295.65	
      yd(4) = 328.56	
      yd(5) = 366.67	
      yd(6) = 245.77	
      yd(7) = 274.88	
      yd(8) = 309.04	
      yd(9) = 349.57	
      yd(10) = 397.34	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 2.040e+3   
      p1noma(2) = 2.801e+3  
      p1noma(3) = 3.766e+3  
      p1noma(4) = 4.956e+3  
      p1noma(5) = 6.386e+3  
      p1noma(6) = 1.635e+3  
      p1noma(7) = 2.336e+3    
      p1noma(8) = 3.249e+3    
      p1noma(9) = 4.394e+3    
      p1noma(10) = 5.787e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 550.6
      fsa(2) = 548.1
      fsa(3) = 557.6
      fsa(4) = 569.1
      fsa(5) = 583.3
      fsa(6) = 428.9
      fsa(7) = 426.7
      fsa(8) = 425.7
      fsa(9) = 426.1
      fsa(10) = 429.9	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s) 
c
      ispa(1) = 1925.4
      ispa(2) = 1917.1
      ispa(3) = 1949.8
      ispa(4) = 1990.0
      ispa(5) = 2040.0
      ispa(6) = 1499.7	  
      ispa(7) = 1492.2
      ispa(8) = 1487.7
      ispa(9) = 1489.7
      ispa(10) = 1503.0	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      phi(1) = 1.00000
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 2.959
      Ma(2) = 2.914
      Ma(3) = 2.916
      Ma(4) = 2.918
      Ma(5) = 2.919
      Ma(6) = 2.995
      Ma(7) = 2.966
      Ma(8) = 2.934
      Ma(9) = 2.900
      Ma(10) = 2.864	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 9802.
      pa(2) = 12141.
      pa(3) = 14254.
      pa(4) = 16270.
      pa(5) = 18067.
      pa(6) = 9796.
      pa(7) = 12281.
      pa(8) = 14903.
      pa(9) = 17456.
      pa(10) = 19786.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1427.9
      ta(2) = 1454.8
      ta(3) = 1452.5
      ta(4) = 1451.3
      ta(5) = 1451.7
      ta(6) = 1599.3
      ta(7) = 1613.4
      ta(8) = 1630.9
      ta(9) = 1650.9
      ta(10) = 1674.0	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.2826
      ga(2) = 1.2813
      ga(3) = 1.2814
      ga(4) = 1.2814
      ga(5) = 1.2814
      ga(6) = 1.2745	
      ga(7) = 1.2739
      ga(8) = 1.2732
      ga(9) = 1.2724
      ga(10) = 1.2715	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.44
      ra(2) = 350.42
      ra(3) = 350.43
      ra(4) = 350.44
      ra(5) = 350.42
      ra(6) = 350.47
      ra(7) = 350.46
      ra(8) = 350.50
      ra(9) = 350.48
      ra(10) = 350.49	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ra, ni, xi, yi, re, iwk, wk )
c
      return
      end  
c
c
c
      subroutine 
     $   bivar_810(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c date - 5 November 2011
c author - M. Smart
c
      integer, parameter :: ni = 1
      integer, parameter :: nd = 10
      integer, parameter :: niwk = 31 * nd + ni
      integer, parameter :: nwk = 8 * nd
c
      integer i
      integer iwk(niwk)
      integer md
      real wk(nwk)
      real xd(nd)
      real xi(ni)
      real yd(nd)
      real yi(ni)
c
      real fsa(nd),ispa(nd),er(nd),Ma(nd),pa(nd),Ta(nd)
      real ga(nd),ra(nd),p1noma(nd)
      real fs(ni),isp(ni),phi(ni),Mache(ni),prese(ni)
      real Tempe(ni),gammae(ni),re(ni),p1nom(ni)
c
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M8; (6-10)M10
c
      xd(1) = 7.534		
      xd(2) = 7.087		
      xd(3) = 6.642	
      xd(4) = 6.199		
      xd(5) = 5.762		
      xd(6) = 9.298		
      xd(7) = 8.623		
      xd(8) = 7.950		
      xd(9) = 7.289	
      xd(10) = 6.657		  
c
      yd(1) = 245.77	
      yd(2) = 274.88	
      yd(3) = 309.04	
      yd(4) = 349.57	
      yd(5) = 397.34	
      yd(6) = 252.59	
      yd(7) = 291.11	
      yd(8) = 338.67	
      yd(9) = 397.36	
      yd(10) = 468.38	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 1.635e+3   
      p1noma(2) = 2.336e+3  
      p1noma(3) = 3.249e+3  
      p1noma(4) = 4.394e+3  
      p1noma(5) = 5.787e+3  
      p1noma(6) = 1.146e+3  
      p1noma(7) = 1.767e+3    
      p1noma(8) = 2.612e+3    
      p1noma(9) = 3.704e+3    
      p1noma(10) = 5.054e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 428.9
      fsa(2) = 426.7
      fsa(3) = 425.7
      fsa(4) = 426.1
      fsa(5) = 429.9
      fsa(6) = 267.8
      fsa(7) = 265.8
      fsa(8) = 259.5
      fsa(9) = 252.6
      fsa(10) = 246.2	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 1499.7
      ispa(2) = 1492.2
      ispa(3) = 1487.7
      ispa(4) = 1489.7
      ispa(5) = 1503.0
      ispa(6) =  934.7	  
      ispa(7) =  929.4
      ispa(8) =  907.2
      ispa(9) =  883.3
      ispa(10) =  860.9	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      phi(1) = 1.00000
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 2.995
      Ma(2) = 2.920
      Ma(3) = 2.934
      Ma(4) = 2.900
      Ma(5) = 2.864
      Ma(6) = 3.396
      Ma(7) = 3.357
      Ma(8) = 3.305
      Ma(9) = 3.241
      Ma(10) = 3.168	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 9796.
      pa(2) = 12281.
      pa(3) = 14903.
      pa(4) = 17456.
      pa(5) = 19786.
      pa(6) = 8232.
      pa(7) = 10923.
      pa(8) = 13681.
      pa(9) = 16315.
      pa(10) = 18682.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1599.3
      ta(2) = 1613.4
      ta(3) = 1630.9
      ta(4) = 1650.9
      ta(5) = 1674.0
      ta(6) = 1721.6
      ta(7) = 1747.2
      ta(8) = 1778.7
      ta(9) = 1819.4
      ta(10) = 1870.3	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.2745
      ga(2) = 1.2739
      ga(3) = 1.2732
      ga(4) = 1.2742
      ga(5) = 1.2715
      ga(6) = 1.2697	
      ga(7) = 1.2688
      ga(8) = 1.2677
      ga(9) = 1.2663
      ga(10) = 1.2647	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.47
      ra(2) = 350.46
      ra(3) = 350.50
      ra(4) = 350.48
      ra(5) = 350.49
      ra(6) = 350.63
      ra(7) = 350.55
      ra(8) = 350.58
      ra(9) = 350.60
      ra(10) = 350.66	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ra, ni, xi, yi, re, iwk, wk )
c
      return
      end  
c
c
c
      subroutine 
     $ bivar_1012(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c date - 5 November 2011
c author - M. Smart
c
      integer, parameter :: ni = 1
      integer, parameter :: nd = 10
      integer, parameter :: niwk = 31 * nd + ni
      integer, parameter :: nwk = 8 * nd
c
      integer i
      integer iwk(niwk)
      integer md
      real wk(nwk)
      real xd(nd)
      real xi(ni)
      real yd(nd)
      real yi(ni)
c
      real fsa(nd),ispa(nd),er(nd),Ma(nd),pa(nd),Ta(nd)
      real ga(nd),ra(nd),p1noma(nd)
      real fs(ni),isp(ni),phi(ni),Mache(ni),prese(ni)
      real Tempe(ni),gammae(ni),re(ni),p1nom(ni)
c
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M10; (6-10)M12
c
      xd(1) = 9.298		
      xd(2) = 8.623		
      xd(3) = 7.950	
      xd(4) = 7.289		
      xd(5) = 6.657		
      xd(6) = 11.011		
      xd(7) = 10.057		
      xd(8) = 9.111		
      xd(9) = 8.204	
      xd(10) = 7.367		  
c
      yd(1) = 252.59	
      yd(2) = 291.11	
      yd(3) = 338.67	
      yd(4) = 397.36	
      yd(5) = 468.38
      yd(6) = 259.66	  
      yd(7) = 308.84	
      yd(8) = 372.43	
      yd(9) = 453.33	
      yd(10) = 552.98	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 1.146e+3   
      p1noma(2) = 1.767e+3  
      p1noma(3) = 2.612e+3  
      p1noma(4) = 3.704e+3  
      p1noma(5) = 5.054e+3  
      p1noma(6) = 0.870e+3  
      p1noma(7) = 1.441e+3    
      p1noma(8) = 2.247e+3    
      p1noma(9) = 3.310e+3    
      p1noma(10) = 4.638e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 267.8
      fsa(2) = 265.8
      fsa(3) = 259.5
      fsa(4) = 252.6
      fsa(5) = 246.2
      fsa(6) = 98.1
      fsa(7) = 108.9
      fsa(8) = 127.6
      fsa(9) = 114.5
      fsa(10) = 99.0	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 934.7	  
      ispa(2) = 929.4
      ispa(3) = 907.2
      ispa(4) = 883.3
      ispa(5) = 860.9	
      ispa(6) = 343.0
      ispa(7) = 380.7
      ispa(8) = 445.9
      ispa(9) = 400.3
      ispa(10) = 346.2	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      phi(1) = 1.00000
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 3.396
      Ma(2) = 3.357
      Ma(3) = 3.305
      Ma(4) = 3.241
      Ma(5) = 3.168
      Ma(6) = 3.817
      Ma(7) = 3.766
      Ma(8) = 3.705
      Ma(9) = 3.620
      Ma(10) = 3.518	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 8232.
      pa(2) = 10923.
      pa(3) = 13681.
      pa(4) = 16315.
      pa(5) = 18682.
      pa(6) = 6659.
      pa(7) = 9466.
      pa(8) = 12585.
      pa(9) = 15094.
      pa(10) = 17272.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1721.6
      ta(2) = 1747.2
      ta(3) = 1778.7
      ta(4) = 1819.4
      ta(5) = 1870.3
      ta(6) = 1798.5
      ta(7) = 1849.7
      ta(8) = 1919.8
      ta(9) = 1977.7
      ta(10) = 2054.1	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.2697
      ga(2) = 1.2688
      ga(3) = 1.2677
      ga(4) = 1.2663
      ga(5) = 1.2647
      ga(6) = 1.2670	
      ga(7) = 1.2654
      ga(8) = 1.2633
      ga(9) = 1.2617
      ga(10) = 1.2599	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.63
      ra(2) = 350.55
      ra(3) = 350.58
      ra(4) = 350.60
      ra(5) = 350.66
      ra(6) = 350.63
      ra(7) = 350.68
      ra(8) = 350.83
      ra(9) = 350.94
      ra(10) = 351.20	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ra, ni, xi, yi, re, iwk, wk )
c
      return
      end  
c
c
c
      subroutine 
     $  bivar_1214(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
c
c***********************************************************************
c
c  Fortran90 Program to interpolate data from an irregular set of data.
c  Uses routine taken from:
c  http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c
c date - 5 November 2011
c author - M. Smart
c
      integer, parameter :: ni = 1
      integer, parameter :: nd = 10
      integer, parameter :: niwk = 31 * nd + ni
      integer, parameter :: nwk = 8 * nd
c
      integer i
      integer iwk(niwk)
      integer md
      real wk(nwk)
      real xd(nd)
      real xi(ni)
      real yd(nd)
      real yi(ni)
c
      real fsa(nd),ispa(nd),er(nd),Ma(nd),pa(nd),Ta(nd)
      real ga(nd),ra(nd),p1noma(nd)
      real fs(ni),isp(ni),phi(ni),Mache(ni),prese(ni)
      real Tempe(ni),gammae(ni),re(ni),p1nom(ni)
c
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M12; (6-10)M14
c
      xd(1) = 11.011	
      xd(2) = 10.057		
      xd(3) = 9.111	
      xd(4) = 8.204		
      xd(5) = 7.367		
      xd(6) = 12.673		
      xd(7) = 11.390		
      xd(8) = 10.133		
      xd(9) = 8.965	
      xd(10) = 7.929		  
c
      yd(1) = 259.66	
      yd(2) = 308.84	
      yd(3) = 372.43	
      yd(4) = 453.33	
      yd(5) = 552.98	
      yd(6) = 267.00	
      yd(7) = 328.22	
      yd(8) = 410.63	
      yd(9) = 517.92	
      yd(10) = 651.62	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 0.870e+3   
      p1noma(2) = 1.441e+3  
      p1noma(3) = 2.247e+3  
      p1noma(4) = 3.310e+3  
      p1noma(5) = 4.638e+3  
      p1noma(6) = 0.698e+3  
      p1noma(7) = 1.236e+3    
      p1noma(8) = 2.017e+3    
      p1noma(9) = 3.063e+3    
      p1noma(10) = 4.379e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 98.1
      fsa(2) = 108.9
      fsa(3) = 127.6
      fsa(4) = 114.5
      fsa(5) = 99.0
      fsa(6) = -48.3
      fsa(7) = -41.4
      fsa(8) = -26.0
      fsa(9) = -12.0
      fsa(10) = -27.1	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 343.0
      ispa(2) = 380.7
      ispa(3) = 445.9
      ispa(4) = 400.3
      ispa(5) = 346.2
      ispa(6) = -168.8	  
      ispa(7) = -144.6
      ispa(8) =  -90.4
      ispa(9) =  -41.9
      ispa(10) = -90.6	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      phi(1) = 1.00000
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 3.871
      Ma(2) = 3.766
      Ma(3) = 3.705
      Ma(4) = 3.620
      Ma(5) = 3.518
      Ma(6) = 4.252
      Ma(7) = 4.182
      Ma(8) = 4.095
      Ma(9) = 4.001
      Ma(10) = 3.848	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 6659.
      pa(2) = 9466.
      pa(3) = 12585.
      pa(4) = 15094.
      pa(5) = 17272.
      pa(6) = 5582.
      pa(7) = 8259.
      pa(8) = 11206.
      pa(9) = 14001.
      pa(10) = 16233.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1798.5
      ta(2) = 1849.7
      ta(3) = 1919.8
      ta(4) = 1977.7
      ta(5) = 2054.1
      ta(6) = 1875.4
      ta(7) = 1939.3
      ta(8) = 2031.4
      ta(9) = 2139.1
      ta(10) = 2257.2	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.2670
      ga(2) = 1.2654
      ga(3) = 1.2633
      ga(4) = 1.2617
      ga(5) = 1.2599
      ga(6) = 1.2647	
      ga(7) = 1.2628
      ga(8) = 1.2604
      ga(9) = 1.2583
      ga(10) = 1.2557	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.63
      ra(2) = 350.68
      ra(3) = 350.83
      ra(4) = 350.94
      ra(5) = 351.20
      ra(6) = 350.80
      ra(7) = 350.92
      ra(8) = 351.43
      ra(9) = 351.84
      ra(10) = 355.88	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ra, ni, xi, yi, re, iwk, wk )
c
      return
      end  
c	
c
c
      subroutine inlet_restm12(M1,ptrati,mc,M2,prati,Trati)
c
      real M1,M2,mc
c
      a1 = 0.5730261
      a2 = 1.65047e-2
      a3 = 8.889289e-3
      a4 = -6.0995262e-4
c
      b1 = -0.25501276
      b2 = 0.40591054
      b3 = 4.5992973e-2
      b4 = -3.3262341e-3
c
      c1 = 8.952760
      c2 = 2.5856736
      c3 = -0.2934742
      c4 = 1.9807402e-2
c
      d1 = 2.0864267
      d2 = -0.67858372
      d3 = 9.780998e-2
      d4 = -4.456572e-3
c
      e1 = 0.6897386
      e2 = 0.7186836
      e3 = -9.1530813e-2
      e4 = 4.659584e-3
c
      mc     = a1 + a2*M1 + a3*M1**2 + a4*M1**3
      M2     = b1 + b2*M1 + b3*M1**2 + b4*M1**3
      prati  = c1 + c2*M1 + c3*M1**2 + c4*M1**3
      ptrati = d1 + d2*M1 + d3*M1**2 + d4*M1**3
      trati  = e1 + e2*M1 + e3*M1**2 + e4*M1**3
c      write(6,*)'mc,M2,prat,ptrat,trat:',
c     $           mc,M2,prati,ptrati,trati
c
      return
      end
c
c