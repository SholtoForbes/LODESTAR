
      program prop

	  real M1,T1,p1,isp,fs,phi,Mache,prese,tempe,gammae,Re

      open(890, file = 'input')
      read(890,*)M1, T1, p1, wcap









      call CRESTM10(M1,T1,p1,isp,fs,phi,Mache,prese,tempe,gammae,Re)

      open(788, file = 'output')
      write(788,*)isp



      end
