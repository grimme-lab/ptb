cccccccccccccccccccccccccccccccccccccccccccccc
c general average for positive values a and b
c e=1  : arithmetic
c e=2  : closer to max value
c e=0.5: closer to min value (like geom. mean)
cccccccccccccccccccccccccccccccccccccccccccccc

      real*8 function gav(e,a,b)
      real*8 e,a,b
      gav = (0.5d0*(a**e+b**e))**(1d0/e)
      end
