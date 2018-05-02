c **********************************
c ***** Program thrust_CREST *******
c **********************************
c
      subroutine thrust_CREST(M0,thetac,alpha,p0,T0,wcap,
     $      fs,f,isp,phi,wf,M1,T1,Mache,prese,tempe,gammae,Re)
c
c  Notes:
c   19 December 2016 - comment out most of the output to screen
c   3 January 2017 - add M1 and T1 as output
c
c
c  Input - M0 - flight Mach number 
c          thetac - forebody cone angle (deg.)
c          alpha - angle of attack (deg.)
c          p0 - atmospheric pressure (Pa)
c          T0 - atmospheric Temperature (K)
c          wcap - capture width of engine (m)
c
c   Output - fs - specific thrust (m/s)
c            f - thrust (N)
c            isp - specific impulse (s)
c            phi - equivalence ratio
c            wf - fuel mass flow (kg/s)
c            M1 - Mach number entering engine
c            T1 - Temperature entering engine
c            Mache - nozzle exit Mach number
c            prese - nozzle exit pressure (Pa)
c            Tempe - nozzle exit Temperature (K)
c            gammae - nozzle exit ratio of specific heats 
c            Re - nozzle exit gas constant (J/kgK)
c
c
c   Date 30 September 2016
c   Author - M. Smart
c
c***********************************************************************************
c  
c 
      real M0,M1,mc,isp,Mrat,Mache
c
c      open(1,file='thrust.in',status='old')
c      read(1,*) M0,thetac,alpha,p0,T0,wcap	  
c 
c     C-REST engine sized for Dawid's baseline vehicle
c
      wcapstandard = 0.5395
c                     meters
c
      gam0 = 1.4000000
      RR = 287.035
c                   J/kgK
c
c   Calculate nominal inflow conditions to engine based on an equivalent cone
c
c      write(*,*)'enter cone shoot'
      call cone_shoot(M0,thetac,alpha,M1,prat,Trat)
c
      p1 = prat*p0
      T1 = Trat*T0
c
c   Determine mass flow through engine
c
      Acapstandard = 0.3057 + 0.06891*(thetac + alpha)
c                              Linear fit
c
      Acap = (wcap/wcapstandard)**2*Acapstandard
c
      Mrat = M0/M1
      call CREST_inlet(M1,Mrat,mc)
c      write(6,*)'M1,mc:',M1,mc
c
      w2 = 0.9*mc*Acap*p0*M0*sqrt(gam0/RR/T0)*4.0 
c                                	  4 engine modules; average of 90% symmetry plane engine over 4 engines
c      w2 = mc*Acap*p0*M0*sqrt(gam0/RR/T0)*4.0 
c                      (kg/s)
c
      call CRESTM10(M1,T1,p1nom,isp,fs,phi,Mache,prese,tempe,gammae,Re)
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
c      write(6,*)'Mach,cone(deg),alpha(deg):',M0,thetac,alpha
c      write(6,*)'p0(Pa),T0(K),wcap(m):',p0,T0,wcap
c      write(6,*)
c      write(6,*)'M1,p1(Pa),T1(K):',M1,p1,T1
c      write(6,*)	
c  
c      write(6,*)'fs(m/s),f(N),isp(s),phi,wf(kg/s):',fs,f,isp,phi,wf
c      write(6,*)'Me,pe(Pa),Te(K),gammae,Re(J/kgK):',
c     $                        Mache,prese,Tempe,gammae,Re
c	  
      end
c
c
c
      subroutine CREST_inlet(M1,Mrat,mc)
c
      real M1,Mrat,MM,mc
c
      MM = 0.175097746892639*M1**(0.776790959520025)
     $               *(Mrat)**(0.57952831191643) - 0.121263698193072
      a1 = -0.08084
      a2 = 0.9422
      a3 = 0.7429
      a4 = -0.6744
c
      mc = a1 + a2*MM + a3*MM**2 + a4*MM**3
c
      return
      end	  
c
c
c
      subroutine CRESTM10(xi,yi,p1nom,isp,fs,phi,
     $                   Mache,prese,tempe,gammae,Re)
c
c  fortran  subroutine to calculate the internal performance of the C-RESTM10 
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
c data - 30 September 2016
c author - M. Smart
c
c notes:  
c
      integer, parameter :: ni = 1
      real xi(ni),yi(ni),e1(ni),e2(ni),e3(ni),e4(ni),e5(ni)
      real e6(ni),e7(ni),e8(ni),e9(ni)
      real isp,Mache
c
      aL =  115.6
      bL =  50.42
      cL =  -6.783
      dL = 0.3881
      tlow  = aL + bL*xi(1) + cL*xi(1)**2 + dL*xi(1)**3
c
      ah = -842.7
      bh =  698.4
      ch = -148.5
      dh =  11.5
      thigh  = ah + bh*xi(1) + ch*xi(1)**2 + dh*xi(1)**3	  
c
      a5 =  1594.
      b5 =   -597.2
      c5 =  88.99
      d5 =  -4.892
      tM5 = a5 + b5*xi(1) + c5*xi(1)**2 + d5*xi(1)**3  
c
      a6 =  2125.
      b6 =   -770.6
      c6 =  110.1
      d6 =  -5.758
      tM6 = a6 + b6*xi(1) + c6*xi(1)**2 + d6*xi(1)**3
c
      a7 =  2273.
      b7 = -704.9
      c7 =  84.77
      d7 = -3.695
      tM7 = a7 + b7*xi(1) + c7*xi(1)**2 + d7*xi(1)**3
c
      a8 =  2593.
      b8 = -742.1
      c8 =  81.76
      d8 = -3.25
      tM8 = a8 + b8*xi(1) + c8*xi(1)**2 + d8*xi(1)**3
c
      a9 =  3046.
      b9 = -836.3
      c9 =  88.19
      d9 = -3.351
      tM9 = a9 + b9*xi(1) + c9*xi(1)**2 + d9*xi(1)**3
c
      a10 =  3284.
      b10 = -823.5
      c10 =  78.75
      d10 = -2.699
      tM10 = a10 + b10*xi(1) + c10*xi(1)**2 + d10*xi(1)**3
c
c      write(6,*) xi,yi
c      write(6,*) tM5,tM6,tM7,tM8,tM9,tM10	  
c
      if(yi(1) .lt. tM5) then	  
        write(6,*)'M1=',xi(1),' below range of propulsion deck'	
        write(6,*)'T1=',yi(1),' tM5=',tM5		
        stop		
      elseif(yi(1) .lt. tM6) then
c          write(6,*)'M1=',xi(1),'Between M5 and M6 lines'
          call bivar_56(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)	
      elseif(yi(1) .lt. tM7) then
c          write(6,*)'M1=',xi(1),'Between M6 and M7 lines'	  
          call bivar_67(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tM8) then
c          write(6,*)'M1=',xi(1),'Between M7 and M8 lines'	  
          call bivar_78(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tM9) then
c          write(6,*)'M1=',xi(1),'Between M8 and M9 lines'	  
          call bivar_89(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)
      elseif(yi(1) .lt. tM10) then
c          write(6,*)'M1=',xi(1),'Between M9 and M10 lines'	  
          call bivar_910(xi,yi,e1,e2,e3,e4,e5,e6,e7,e8,e9)	
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
      return
      end
c
c
c
      subroutine 
     $   bivar_56(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
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
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M5; (6-10)M6
c
      xd(1) = 4.71478		
      xd(2) = 4.55187		
      xd(3) = 4.38037	
      xd(4) = 4.20261		
      xd(5) = 4.02037		
      xd(6) = 5.58210		
      xd(7) = 5.35601		
      xd(8) = 5.11583		
      xd(9) = 4.86831	
      xd(10) = 4.61730		  
c
      yd(1) = 243.269	
      yd(2) = 257.548	
      yd(3) = 273.858	
      yd(4) = 292.295	
      yd(5) = 312.993	
      yd(6) = 252.834	
      yd(7) = 271.654	
      yd(8) = 293.573	
      yd(9) = 318.852	
      yd(10) = 347.696	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 4.00943e+3   
      p1noma(2) = 4.87778e+3  
      p1noma(3) = 5.96673e+3  
      p1noma(4) = 7.27881e+3  
      p1noma(5) = 8.81365e+3  
      p1noma(6) = 3.06522e+3  
      p1noma(7) = 3.89973e+3    
      p1noma(8) = 4.95984e+3    
      p1noma(9) = 6.24781e+3    
      p1noma(10) = 7.76264e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 515.13
      fsa(2) = 510.07
      fsa(3) = 486.72
      fsa(4) = 472.85
      fsa(5) = 448.06
      fsa(6) = 662.40
      fsa(7) = 666.92
      fsa(8) = 671.89
      fsa(9) = 651.80
      fsa(10) = 586.01	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 2776.13
      ispa(2) = 2791.84
      ispa(3) = 2795.04
      ispa(4) = 2807.42
      ispa(5) = 2753.61
      ispa(6) = 2320.38	  
      ispa(7) = 2336.20
      ispa(8) = 2353.63
      ispa(9) = 2283.26
      ispa(10) = 2280.88	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      er(1) = 0.65
      er(2) = 0.64
      er(3) = 0.61
      er(4) = 0.59
      er(5) = 0.57
      er(6) = 1.0
      er(7) = 1.0
      er(8) = 1.0
      er(9) = 1.0
      er(10) = 0.90 
cc  to here
cc
      md = 3
c 
      call idbvip ( md, nd, xd, yd, er, ni, xi, yi, phi, iwk, wk )
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 2.687
      Ma(2) = 2.688
      Ma(3) = 2.691
      Ma(4) = 2.693
      Ma(5) = 2.697
      Ma(6) = 2.694
      Ma(7) = 2.695
      Ma(8) = 2.696
      Ma(9) = 2.665
      Ma(10) = 2.676	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 11457.
      pa(2) = 12842.
      pa(3) = 14128.
      pa(4) = 15601.
      pa(5) = 16914.
      pa(6) = 13408.
      pa(7) = 15480.
      pa(8) = 17781.
      pa(9) = 20469.
      pa(10) = 22114.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1066.8
      ta(2) = 1058.7
      ta(3) = 1034.4
      ta(4) = 1017.9
      ta(5) = 992.9
      ta(6) = 1477.6
      ta(7) = 1476.0
      ta(8) = 1475.0
      ta(9) = 1465.6
      ta(10) = 1397.0	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.314
      ga(2) = 1.315
      ga(3) = 1.317
      ga(4) = 1.319
      ga(5) = 1.322
      ga(6) = 1.280	
      ga(7) = 1.280
      ga(8) = 1.280
      ga(9) = 1.282
      ga(10) = 1.287	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit 
c
      ra(1) = 330.3
      ra(2) = 329.7
      ra(3) = 327.9
      ra(4) = 326.6
      ra(5) = 324.8
      ra(6) = 350.5
      ra(7) = 350.4
      ra(8) = 350.4
      ra(9) = 351.6
      ra(10) = 345.6	  
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
     $  bivar_67(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
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
      xd(1) = 5.58210		
      xd(2) = 5.35601		
      xd(3) = 5.11583		
      xd(4) = 4.86831	
      xd(5) = 4.61730		
      xd(6) = 6.42687		
      xd(7) = 6.11801		
      xd(8) = 5.79580		
      xd(9) = 5.46701	
      xd(10) = 5.13886		  
c
      yd(1) = 252.834	
      yd(2) = 271.654	
      yd(3) = 293.573	
      yd(4) = 318.852	
      yd(5) = 347.696	
      yd(6) = 262.626	
      yd(7) = 286.608	
      yd(8) = 315.118	
      yd(9) = 348.565	
      yd(10) = 387.189	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 3.06522e+3  
      p1noma(2) = 3.89973e+3    
      p1noma(3) = 4.95984e+3    
      p1noma(4) = 6.24781e+3    
      p1noma(5) = 7.76264e+3 
      p1noma(6) = 2.48755e+3  
      p1noma(7) = 3.30055e+3    
      p1noma(8) = 4.34341e+3    
      p1noma(9) = 5.61782e+3    
      p1noma(10) = 7.12199e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 662.40
      fsa(2) = 666.92
      fsa(3) = 671.89
      fsa(4) = 651.80
      fsa(5) = 586.01
      fsa(6) = 631.32
      fsa(7) = 633.08
      fsa(8) = 504.20
      fsa(9) = 509.23
      fsa(10) = 519.34	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s) 
c
      ispa(1) = 2320.38	  
      ispa(2) = 2336.20
      ispa(3) = 2353.63
      ispa(4) = 2283.26
      ispa(5) = 2280.88
      ispa(6) = 2211.50	  
      ispa(7) = 2217.65
      ispa(8) = 1766.21
      ispa(9) = 1783.83
      ispa(10) = 1819.22	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ispa, ni, xi, yi, isp, iwk, wk )
c
c   equivalence ratio
c
      er(1) = 1.0
      er(2) = 1.0
      er(3) = 1.0
      er(4) = 1.0
      er(5) = 0.9
      er(6) = 1.0
      er(7) = 1.0
      er(8) = 1.0
      er(9) = 1.0
      er(10) = 1.0 
cc  to here
cc
      md = 3
c 
      call idbvip ( md, nd, xd, yd, er, ni, xi, yi, phi, iwk, wk )
c
c  Mach number at Internal nozzle exit  
c
      Ma(1) = 2.694
      Ma(2) = 2.695
      Ma(3) = 2.696
      Ma(4) = 2.665
      Ma(5) = 2.676
      Ma(6) = 2.717
      Ma(7) = 2.707
      Ma(8) = 2.697
      Ma(9) = 2.689
      Ma(10) = 2.692	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 13408.
      pa(2) = 15480.
      pa(3) = 17781.
      pa(4) = 20469.
      pa(5) = 22114.
      pa(6) = 13374.
      pa(7) = 15938.
      pa(8) = 18936.
      pa(9) = 22016.
      pa(10) = 24771.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1477.6
      ta(2) = 1476.0
      ta(3) = 1475.0
      ta(4) = 1465.6
      ta(5) = 1397.0
      ta(6) = 1629.9
      ta(7) = 1635.4
      ta(8) = 1641.3
      ta(9) = 1645.3
      ta(10) = 1643.5	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.280	
      ga(2) = 1.280
      ga(3) = 1.280
      ga(4) = 1.282
      ga(5) = 1.287
      ga(6) = 1.273	
      ga(7) = 1.273
      ga(8) = 1.273
      ga(9) = 1.273
      ga(10) = 1.273	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.5
      ra(2) = 350.4
      ra(3) = 350.4
      ra(4) = 351.6
      ra(5) = 345.6
      ra(6) = 350.5
      ra(7) = 350.5
      ra(8) = 350.5
      ra(9) = 350.5
      ra(10) = 350.5	  
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
     $   bivar_78(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
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
      xd(1) = 6.42687		
      xd(2) = 6.11801		
      xd(3) = 5.79580		
      xd(4) = 5.46701	
      xd(5) = 5.13886		  	  
      xd(6) = 7.23904		
      xd(7) = 6.83691		
      xd(8) = 6.42035		
      xd(9) = 6.00128	
      xd(10) = 5.59104			  
c
      yd(1) = 262.626	
      yd(2) = 286.608	
      yd(3) = 315.118	
      yd(4) = 348.565	
      yd(5) = 387.189		
      yd(6) = 272.616	
      yd(7) = 302.438	
      yd(8) = 338.574	
      yd(9) = 381.564	
      yd(10) = 431.587
c	  
c  to here
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 2.48755e+3  
      p1noma(2) = 3.30055e+3    
      p1noma(3) = 4.34341e+3    
      p1noma(4) = 5.61782e+3    
      p1noma(5) = 7.12199e+3   
      p1noma(6) = 2.10614e+3  
      p1noma(7) = 2.90438e+3    
      p1noma(8) = 3.93585e+3    
      p1noma(9) = 5.20152e+3    
      p1noma(10) = 6.69883e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 631.32
      fsa(2) = 633.08
      fsa(3) = 504.20
      fsa(4) = 509.23
      fsa(5) = 519.34
      fsa(6) = 391.61
      fsa(7) = 391.39
      fsa(8) = 391.21
      fsa(9) = 392.94
      fsa(10) = 412.68	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 2211.50	  
      ispa(2) = 2217.65
      ispa(3) = 1766.21
      ispa(4) = 1783.83
      ispa(5) = 1819.22
      ispa(6) =  1371.82	  
      ispa(7) =  1371.02
      ispa(8) =  1370.40
      ispa(9) =  1376.45
      ispa(10) =  1445.61	  
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
      Ma(1) = 2.717
      Ma(2) = 2.707
      Ma(3) = 2.697
      Ma(4) = 2.689
      Ma(5) = 2.692
      Ma(6) = 2.864
      Ma(7) = 2.845
      Ma(8) = 2.821
      Ma(9) = 2.793
      Ma(10) = 2.750	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 13374.
      pa(2) = 15938.
      pa(3) = 18936.
      pa(4) = 22016.
      pa(5) = 24771.
      pa(6) = 12461.
      pa(7) = 15631.
      pa(8) = 18981.
      pa(9) = 22334.
      pa(10) = 25420.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1629.9
      ta(2) = 1635.4
      ta(3) = 1641.3
      ta(4) = 1645.3
      ta(5) = 1643.5
      ta(6) = 1732.9
      ta(7) = 1744.7
      ta(8) = 1758.8
      ta(9) = 1773.4
      ta(10) = 1772.1	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.273	
      ga(2) = 1.273
      ga(3) = 1.273
      ga(4) = 1.273
      ga(5) = 1.273
      ga(6) = 1.269	
      ga(7) = 1.269
      ga(8) = 1.268
      ga(9) = 1.268
      ga(10) = 1.268	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.5
      ra(2) = 350.5
      ra(3) = 350.5
      ra(4) = 350.5
      ra(5) = 350.5
      ra(6) = 350.5
      ra(7) = 350.5
      ra(8) = 350.5
      ra(9) = 350.5
      ra(10) = 350.5	  
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
     $ bivar_89(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
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
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M8; (6-10)M9
c
      xd(1) = 7.23904		
      xd(2) = 6.83691		
      xd(3) = 6.42035		
      xd(4) = 6.00128	
      xd(5) = 5.59104			
      xd(6) = 8.02107		
      xd(7) = 7.51232		
      xd(8) = 6.99075		
      xd(9) = 6.47522	
      xd(10) = 5.98110		  
c
      yd(1) = 272.616	
      yd(2) = 302.438	
      yd(3) = 338.574	
      yd(4) = 381.564	
      yd(5) = 431.587
      yd(6) = 283.411	  
      yd(7) = 319.866	
      yd(8) = 364.782	
      yd(9) = 418.745	
      yd(10) = 481.954	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 2.10614e+3  
      p1noma(2) = 2.90438e+3    
      p1noma(3) = 3.93585e+3    
      p1noma(4) = 5.20152e+3    
      p1noma(5) = 6.69883e+3 
      p1noma(6) = 1.84208e+3  
      p1noma(7) = 2.63046e+3    
      p1noma(8) = 3.65487e+3    
      p1noma(9) = 4.91533e+3    
      p1noma(10) = 6.40925e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 391.61
      fsa(2) = 391.39
      fsa(3) = 391.21
      fsa(4) = 392.94
      fsa(5) = 412.68
      fsa(6) = 318.14
      fsa(7) = 308.87
      fsa(8) = 307.38
      fsa(9) = 300.29
      fsa(10) = 292.86	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) =  1371.82	  
      ispa(2) =  1371.02
      ispa(3) =  1370.40
      ispa(4) =  1376.45
      ispa(5) =  1445.61	
      ispa(6) = 1115.38
      ispa(7) = 1081.98
      ispa(8) = 1076.73
      ispa(9) = 1051.93
      ispa(10) = 1025.87	  
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
      Ma(1) = 2.864
      Ma(2) = 2.845
      Ma(3) = 2.821
      Ma(4) = 2.793
      Ma(5) = 2.750
      Ma(6) = 3.073
      Ma(7) = 3.046
      Ma(8) = 3.015
      Ma(9) = 2.978
      Ma(10) = 2.932	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 12461.
      pa(2) = 15631.
      pa(3) = 18981.
      pa(4) = 22334.
      pa(5) = 25420.
      pa(6) = 11733.
      pa(7) = 14900.
      pa(8) = 18270.
      pa(9) = 20412.
      pa(10) = 20212.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1732.9
      ta(2) = 1744.7
      ta(3) = 1758.8
      ta(4) = 1773.4
      ta(5) = 1772.1
      ta(6) = 1803.0
      ta(7) = 1820.0
      ta(8) = 1839.8
      ta(9) = 1864.3
      ta(10) = 1894.3	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.269	
      ga(2) = 1.269
      ga(3) = 1.268
      ga(4) = 1.268
      ga(5) = 1.268
      ga(6) = 1.270	
      ga(7) = 1.266
      ga(8) = 1.266
      ga(9) = 1.265
      ga(10) = 1.264	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.5
      ra(2) = 350.5
      ra(3) = 350.5
      ra(4) = 350.5
      ra(5) = 350.5
      ra(6) = 350.6
      ra(7) = 350.6
      ra(8) = 350.6
      ra(9) = 350.6
      ra(10) = 350.7	  
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
     $  bivar_910(xi,yi,p1nom,fs,isp,phi,Mache,prese,tempe,gammae,Re)
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
c     Independent variables  M1 --- xd;  T1 -- yd (1-5)M9; (6-10)M10
c
      xd(1) = 8.02107		
      xd(2) = 7.51232		
      xd(3) = 6.99075		
      xd(4) = 6.47522	
      xd(5) = 5.98110		
      xd(6) = 8.77245		
      xd(7) = 8.14441		
      xd(8) = 7.50927		
      xd(9) = 6.89391	
      xd(10) = 6.31673		  
c
      yd(1) = 283.411	  
      yd(2) = 319.866	
      yd(3) = 364.782	
      yd(4) = 418.745	
      yd(5) = 481.954		
      yd(6) = 298.002	
      yd(7) = 342.387	
      yd(8) = 397.893	
      yd(9) = 464.970	
      yd(10) = 543.928	
c
c  Data to be interpolated from
c
c  nominal P1 (kPa) 
c
      p1noma(1) = 1.84208e+3  
      p1noma(2) = 2.63046e+3    
      p1noma(3) = 3.65487e+3    
      p1noma(4) = 4.91533e+3    
      p1noma(5) = 6.40925e+3 
      p1noma(6) = 1.64961e+3  
      p1noma(7) = 2.43073e+3    
      p1noma(8) = 3.44922e+3    
      p1noma(9) = 4.70566e+3    
      p1noma(10) = 6.19596e+3  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, p1noma, ni, xi, yi, p1nom, iwk, wk )
c
c  specific thrust (m/s) 
c
      fsa(1) = 318.14
      fsa(2) = 308.87
      fsa(3) = 307.38
      fsa(4) = 300.29
      fsa(5) = 292.86
      fsa(6) = 236.61
      fsa(7) = 235.33
      fsa(8) = 228.30
      fsa(9) = 217.02
      fsa(10) = 206.23	  
c
      md = 1
c 
      call idbvip ( md, nd, xd, yd, fsa, ni, xi, yi, fs, iwk, wk )
c
c  specific Impulse (s)  
c
      ispa(1) = 1115.38
      ispa(2) = 1081.98
      ispa(3) = 1076.73
      ispa(4) = 1051.93
      ispa(5) = 1025.87	
      ispa(6) = 828.85	  
      ispa(7) = 824.35
      ispa(8) = 799.75
      ispa(9) = 760.22
      ispa(10) = 722.41	  
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
      Ma(1) = 3.073
      Ma(2) = 3.046
      Ma(3) = 3.015
      Ma(4) = 2.978
      Ma(5) = 2.932
      Ma(6) = 3.303
      Ma(7) = 3.272
      Ma(8) = 3.231
      Ma(9) = 3.179
      Ma(10) = 3.122	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, Ma, ni, xi, yi, Mache, iwk, wk )
c
c  Pressure at Internal nozzle exit  
c
      pa(1) = 11733.
      pa(2) = 14900.
      pa(3) = 18270.
      pa(4) = 20412.
      pa(5) = 20212.
      pa(6) = 10817.
      pa(7) = 13966.
      pa(8) = 16173.
      pa(9) = 17952.
      pa(10) = 19557.	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, pa, ni, xi, yi, prese, iwk, wk )
c
c  Temperature at Internal nozzle exit  
c
      ta(1) = 1803.0
      ta(2) = 1820.0
      ta(3) = 1839.8
      ta(4) = 1864.3
      ta(5) = 1894.3	
      ta(6) = 1875.6
      ta(7) = 1895.7
      ta(8) = 1922.9
      ta(9) = 1958.2
      ta(10) = 1966.5	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ta, ni, xi, yi, tempe, iwk, wk )
c
c  Gamma at Internal nozzle exit  
c
      ga(1) = 1.270	
      ga(2) = 1.266
      ga(3) = 1.266
      ga(4) = 1.265
      ga(5) = 1.264
      ga(6) = 1.265	
      ga(7) = 1.264
      ga(8) = 1.263
      ga(9) = 1.262
      ga(10) = 1.261	  
c
      md = 3
c 
      call idbvip ( md, nd, xd, yd, ga, ni, xi, yi, gammae, iwk, wk )
c
c  Gas Constant at Internal nozzle exit  
c
      ra(1) = 350.6
      ra(2) = 350.6
      ra(3) = 350.6
      ra(4) = 350.6
      ra(5) = 350.7	
      ra(6) = 350.7
      ra(7) = 350.7
      ra(8) = 351.8
      ra(9) = 351.8
      ra(10) = 355.9	  
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
      subroutine cone_shoot(M1,thetac,alpha,Mc,pcop1,TcoT1)
c
c    M1 - flight Mach number
c    thetac - cone half angle (degrees)
c    alpha - angle-of-attack (degrees)
c
c      PROGRAM cone-shoot
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
c      write(6,'(a,f8.2,a,f5.2)')'  M1=',Mach,
c     $             '  Equiv. cone half angle=',hangle
c      write(6,*)'  Mach=',Mach,'  Equiv. cone half angle=',hangle
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
c        write(6,*) 'Cone: Vrat,Mach,pc/p1,Tc/T1',Vrc,Mc,pcop1,TcoT1
c
      return
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
        write(6,*)'roots not bracketed ??'
        write(6,*)'fa,fb:',fa,fb 
        pause
     *'root not bracketed by initial limits ??'
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
c      write(6,*)' Shock angle, Vtheta on cone:',x1*180./pi,vthetac
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
c        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        if(xnew.eq.x) write(*,*)'stepsize underflow in rkqs'		
        goto 1
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
