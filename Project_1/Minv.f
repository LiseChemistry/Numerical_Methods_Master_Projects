!     Exemple d'utilisation de Minv
      Implicit none
      Integer :: i,j,k,nX
      Integer,dimension(:),allocatable :: L
      Integer,parameter :: dimA=20
      Real*8 :: dfloat,tol
      Real*8,dimension(:),allocatable :: Wk
      Real*8,dimension(dimA,dimA) :: A,Am1
!
!     Construction de la matrice A_{ij} telle que
!          A_{ii}=1
!          A_{ij}=1/(i-j)^2 si i/=j
      nX=8
      do i=1,nX
         A(i,i)=1d0
         do j=1,i-1
            A(i,j)=0.5d0/dfloat(i-j)**2
            A(j,i)=A(i,j)
         enddo
      enddo
      Am1(1:nX,1:nX)=A(1:nX,1:nX)
      Write (6,"('A[i,j]')")
      Write (6,"(8f10.6)") ((A(i,j),i=1,nX),j=1,nX)
!
!     Inversion de la matrice
      Allocate(L(2*nX),Wk(2*nX))
      tol=0d0
      call Minv(Am1,nX,dimA,L,Wk,tol)
      Write (6,"(/'det[A]=',d12.5)") Wk(1)
      Write (6,"(/'A[i,j]^{-1}')")
      Write (6,"(8f10.6)") ((Am1(i,j),i=1,nX),j=1,nX)
!
!     Vérification A.A^{-1}=1
      Write (6,"(/'A.A^{-1}')")
      do i=1,nX
         Wk(1:nX)=0d0
         do j=1,nX
            do k=1,nX
               Wk(j)=Wk(j)+A(i,k)*Am1(k,j)
            enddo
         enddo
         Write (6,"(8f10.6)") Wk(1:nX)
      enddo
      End
!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine Minv(A,nA,dimA,L,Wk,tol)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Matrix inversion and solution of linear systems
!     A(dimA,1) : matrix to be inverted
!     dimA : dimension declared in the calling sequence
!     nA : dimension of the system
!     L    : work array of dimension 2*n
!     Wk   :  "     "   "      "     " .On return, Wk(1) contains the
c            determinant.
c     Tol : tolerence on pivoting (to be set to 0d0 in the call)
!     ===+=========+=========+=========+=========+=========+=========+==
      Implicit Real*8(a-h,o-z)
      Integer dimA
      Dimension A(*),L(*),Wk(*)
c        search for largest element
c
      d=1d0
      nk=-dima
      do 80 k=1,nA
      nk=nk+dima
      l(k)=k
      l(nA+k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,nA
      iz=dima*(j-1)
      do 20 i=k,nA
      ij=iz+i
   10 if(abs(biga)-abs(a(ij))) 15,20,20
   15 biga=a(ij)
      l(k)=i
      l(nA+k)=j
   20 continue
c
c        interchange rows
c
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-dima
      do 30 i=1,nA
      ki=ki+dima
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
c
c        interchange columns
c
   35 i=l(nA+k)
      if(i-k) 45,45,38
   38 jp=dima*(i-1)
      do 40 j=1,nA
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
c
c        divide column by
c        contained in biga)
c
   45 if(abs(biga).le.tol) then
          wk(1)=0.
          return
      endif
      do 55 i=1,nA
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
c
c        reduce matrix
c
      do 65 i=1,nA
      ik=nk+i
      hold=a(ik)
      ij=i-dima
      do 65 j=1,nA
      ij=ij+dima
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue
c
c        divide row by pivot
c
      kj=k-dima
      do 75 j=1,nA
      kj=kj+dima
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
c
c        product of pivots
      d=d*biga
c
c        replace pivot by reciprocal
      a(kk)=1.d0/biga
   80 continue
c
c        final row and column interchange
      k=nA
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=dima*(k-1)
      jr=dima*(i-1)
      do 110 j=1,nA
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=l(nA+k)
      if(j-k) 100,100,125
  125 ki=k-dima
      do 130 i=1,nA
      ki=ki+dima
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
c
  150 Wk(1)=d
      Return
      End
