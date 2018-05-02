      program prop

	  real M1,p1,T1,wcap,fs,f,isp,phi,wf,Mache,prese,tempe,gammae,Re

      M1 = 6
      T1 = 10
      p1 = 50000
      wcap = 1

      call thrust_5_int (M1,p1,T1,wcap,fs,f,isp,phi,wf,Mache,prese,tempe,gammae,Re) ! this needs to have a space, dont know why

      end




!in fixed form things need to be written 7 characters in because thats how they did it on punch cards...

