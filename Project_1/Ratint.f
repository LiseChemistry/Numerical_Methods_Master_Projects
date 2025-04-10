      Subroutine ratint(xa,ya,n,x,y,dy)
c     ============================================================
c     xa : tableau des abcisses de la grille
c     ya : valeurs correspondantes de la fonction
c      n : nombre de points
c      x : abcisse d'interpolation
c      y : valeur interpolée en x
c     dy : estimation de l'erreur d'interpolation
c     ============================================================
      Implicit Real*8(a-h,o-z)
      PARAMETER (NMAX=50,TINY=1.d-25)
      Dimension xa(n),ya(n),c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0d0)stop 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
