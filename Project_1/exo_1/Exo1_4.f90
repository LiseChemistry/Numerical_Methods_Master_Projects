!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine SplSet(F,X,NX,dF1,dFn)
!     ===+=========+=========+=========+=========+=========+=========+==
!         THIS ROUTINE SETS THE SPLINE INTERPOLATION ON THE GRID
!     (X(I),I=1,NX) FOR THE FUNCTION (F(I),I=1,NX).THE F-ARRAY MUST
!     BE DIMENSIONED TO 2*NX .
!     ===+=========+=========+=========+=========+=========+=========+==
      IMPLICIT Real*8 (A-H,O-Z)
      Real*8 :: F(*),X(*)
      Real*8,allocatable,dimension(:) :: HX,RLX,RMUX,B
      DATA UN/1D0/,THREE/3D0/
!
      Allocate(HX(nX),RLX(nX),RMUX(nX),B(nX))
      NX2=NX-2
!     COMPUTE THE HX'S
      DO 10 I=2,NX
   10 HX(I-1)=X(I)-X(I-1)
!     COMPUTE LAMBDA'S & MU'S
      DO 40 I=1,NX2
      RLX(I)=HX(I+1)/(HX(I)+HX(I+1))
   40 RMUX(I)=UN-RLX(I)
!
!     SPLINE-FIT DE P(X)
      F(NX+1)=dF1
      F(2*NX)=dFn
!     CALCUL SECOND MEMBRE
      DO 120 I=1,NX2
  120 B(I)=THREE*RLX(I)/HX(I)*(F(I+1)-F(I))+THREE*RMUX(I)/HX(I+1)*(F(I+2)-F(I+1))
      B(1)=B(1)-RLX(1)*F(NX+1)
      B(NX2)=B(NX2)-RMUX(NX2)*F(2*NX)
      CALL JORDAN(RMUX,RLX,F(NX+2),NX2,B)
      Deallocate(HX,RLX,RMUX,B)
      RETURN
      END
!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine Jordan(mu,lambda,x,n,b)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Résolution du système linéaire
!     ===+=========+=========+=========+=========+=========+=========+==
      IMPLICIT Real*8 (A-H,O-Z)
      Real*8 :: MU(N),LAMBDA(N),X(N),B(N)
      Real*8,allocatable,dimension(:) :: piv
!
!     CALCUL DES PIVOTS
      Allocate(piv(n))
      PIV(1)=2D0
      DO 10 I=2,N
      PIV(I)=2D0-LAMBDA(I)*MU(I-1)/PIV(I-1)
   10 B(I)=B(I)-LAMBDA(I)/PIV(I-1)*B(I-1)
!
      X(N)=B(N)/PIV(N)
      I=N-1
   20 X(I)=(B(I)-X(I+1)*MU(I))/PIV(I)
      I=I-1
      IF(I.GT.0) GOTO 20
      Deallocate(piv)
      RETURN
      END
!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine Splint(F,x,nX,R,F_R,dF_R)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Evaluation de la valeur interpolée F_R en R, et de sa dérivée dF_R
!     ===+=========+=========+=========+=========+=========+=========+==
      Implicit Real*8(A-H,O-Z)
      Dimension F(*),X(*),U(4)
      DATA UN/1D0/,TWO/2D0/,THREE/3D0/
!
!     Linear extrapolation if outside the interpolation range
      if(R.lt.x(1)) then
!     f(r)=f1+(r-r1)*p1=(f1-r1.p1)+r.p1
         F_R=F(1)+(R-x(1))*F(nX+1)  ! =p_1
         goto 999
      Endif
      if(R.gt.x(nX)) Then
!     f(r)=fn+(r-rn)*pn=(fn-rn.pn)+r.pn
         F_R=F(nX)+(R-x(nX))*F(2*nX)  ! =p_N
         goto 999
      Endif
      do IS=2,NX
         if(R.gt.x(iS)) Cycle
         HI=X(IS)-X(IS-1)
         XR=(R-X(IS-1))/HI
         u(1)=XR*XR*(-TWO*XR+THREE)
         u(3)=HI*XR*XR*(XR-UN)
         u(2)=UN-U(1)
         u(4)=HI*XR*((XR-TWO)*XR+UN)
         F_R=U(1)*F(IS)+U(2)*F(IS-1)+U(3)*F(NX+IS)+U(4)*F(NX+IS-1)
!
!     Derivative
         u(1)=-6*(XR-1)*XR/HI
         u(2)=-u(1)
         u(3)=(3*XR-2)*XR
         u(4)=u(3)-2*XR+1
         dF_R=u(1)*F(IS)+u(2)*F(IS-1)+u(3)*F(NX+IS)+u(4)*F(NX+IS-1)
         exit
      enddo
  999 Return
      End Subroutine Splint

program splines_cubiques
  implicit Real*8 (A-H,O-Z)
  real*8 :: D, Rast, alpha, c1, c2, c3
 ! real*8 :: dF1, dFn, x, F_R, dF_R
  integer, parameter :: n = 8
  real*8,dimension(n) :: R = [1.8d0, 2.1d0, 2.4d0, 2.8d0, 3.3d0, 4.0d0, 5.0d0, 7.0d0]
  real*8,dimension(2*n) :: f
  !integer :: i

  D = 3.8886d0
  Rast = 2.2818d0
  alpha = 3.3522498d0
  c1 = 3.6445906d0
  c2 = 3.9281238d0
  c3 = 2.0986689d0
  

  do i = 1, n
    f(i)= -D*(1.0+c1*(R(i)-Rast)+c2*(R(i)-Rast)**2+c3*(R(i)-Rast)**3)*DEXP(-alpha*(R(i)-Rast))
  end do

      dF1=-D*((c1+2*c2*(R(1)-Rast)+3*c3*(R(1)-Rast)**2)*DEXP(-alpha*(R(1)-Rast)) & 
      -alpha*(1+c1*(R(1)-Rast)+c2*(R(1)-Rast)**2+c3*(R(1)-Rast)**3)*DEXP(-alpha*(R(1)-Rast)))
      print*, R(1), dF1
      dFn=-D*((c1+2*c2*(R(n)-Rast)+3*c3*(R(n)-Rast)**2)*DEXP(-alpha*(R(n)-Rast)) & 
      -alpha*(1+c1*(R(n)-Rast)+c2*(R(n)-Rast)**2+c3*(R(n)-Rast)**3)*DEXP(-alpha*(R(n)-Rast)))
      print*, R(n), dFn
!
!     Initialisation (appelée une seule fois)
      call SplSet(f,R,n,dF1,dFn)
!     Interpolation
      do i=1,100
         x=0.1d0*i
         call Splint(f,R,n,x,F_R,dF_R)
         open(unit = 10, file = 'potentiel.dat', status = 'unknown')
         write(10,*) x, F_R
      enddo
 call system ('gnuplot -p potentiel.plt')
 stop
end program splines_cubiques