
      program prop

	  real M1,p1,T1,wcap,fs,f,isp,phi,wf,Mache,prese,tempe,gammae,Re
      
      open(890, file = 'input')
      read(890,*)M1, T1, p1, wcap

    


      




      call thrust_5_int(M1,p1,T1,wcap,fs,f,isp,phi,wf,Mache,
     $                  prese,tempe,gammae,Re)

      open(788, file = 'output')
      write(788,*)f, wf



      end
