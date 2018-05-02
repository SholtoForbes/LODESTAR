
      program prop

	  real M1,thetac,alpha,Mc,pcop1,TcoT1

      open(890, file = 'input')
      read(890,*)M1,thetac,alpha









      call cone_shoot(M1,thetac,alpha,Mc,pcop1,TcoT1)

      open(788, file = 'output')
      write(788,*)Mc,pcop1,TcoT1



      endprogram
