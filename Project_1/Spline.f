      Implicit Real*8 (A-H,O-Z)
      Integer,parameter :: nbX=20
      Real*8,dimension(nbX) :: X
!     Attn : dimension de F est 2 fois le nombre de points
      Real*8,dimension(2*nbX) :: F
!
!     Exemple de Fonction sin(x)**2 à interpoler
      Fonc(t)=sin(t)**2
!     Dérivée analytique
      dFonc(t)=2*sin(t)*cos(t)
!
!     Définition des couples {x_i,F_i} à interpoler
      twopi=2*dacos(-1d0)
      nX=10
      dX=twopi/(nX-1)
!     Grid setup
      x(1)=0d0
      F(1)=Fonc(x(1))
      do i=2,nX
         x(i)=x(i-1)+dX
         F(i)=Fonc(x(i))
      enddo
!
!     Définition des dérivées aux bornes, à modifier pour le potentiel
!     de O2
      dF1=0d0
      dFn=0d0
!
!     Initialisation (appelée une seule fois)
      call SplSet(F,x,nX,dF1,dFn)
!
!     Test d'interpolation
      Write (6,"(13x,'Exact',4x,'Interpolée')")
      do i=1,nX-1
         r=x(i)+dX/2
         call Splint(F,x,nX,r,F_R,dF_R)
         Write (6,"('x=',f6.3)") r
         Write (6,"(5x,'f =',2f12.4)") Fonc(r),F_R
         Write (6,"(5x,'f''=',2f12.4)") dFonc(r),dF_R
      enddo
      stop
      End
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
c
      Allocate(HX(nX),RLX(nX),RMUX(nX),B(nX))
      NX2=NX-2
!     COMPUTE THE HX'S
      DO 10 I=2,NX
   10 HX(I-1)=X(I)-X(I-1)
!     COMPUTE LAMBDA'S & MU'S
      DO 40 I=1,NX2
      RLX(I)=HX(I+1)/(HX(I)+HX(I+1))
   40 RMUX(I)=UN-RLX(I)
C
!     SPLINE-FIT DE P(X)
      F(NX+1)=dF1
      F(2*NX)=dFn
!     CALCUL SECOND MEMBRE
      DO 120 I=1,NX2
  120 B(I)=THREE*RLX(I)/HX(I)*(F(I+1)-F(I))
     & +THREE*RMUX(I)/HX(I+1)*(F(I+2)-F(I+1))
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
C
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
