program main
      implicit none
      integer :: i, m, n
      real(kind=8), dimension(24) :: wa, xa
      real(kind=8) :: ro, r3s, sum
  
      do n = 4, 24, 4
          call laguerreq(n, xa, wa)
          print *, ""
          print *, ""
          print *, ""
          write(*, '(A,I0,A)') "Quadrature à ", n, " points :"
          do i = 1, n
              write(*, '(2E26.17)') xa(i), wa(i)
          end do
  
          ! Vérification \int_0^\infty x^p.e^{-x}dx=p!
          sum = 0.0d0
          r3s = 0.0d0
          do m = 1, n
              sum = sum + wa(m) * (xa(m)**2 - 2*xa(m) + 3)
              ro = xa(m) * 3 / 2
              r3s = r3s + wa(m) * ((2 / sqrt(3.0d0) / 81 * (27 - 18*ro + 2*ro**2))**2) * ro**3
          end do
          write(*, '(A,E26.17,A)') "Vérification :", sum, ", vs 3"
          write(*, '(A,E26.17,A)') "<r>_3s :", r3s * 3 / 2, ", vs. 13.5"
      end do
  end program main
  
  subroutine laguerreq(nx, x_q, w_q)
      implicit none
      integer, intent(inout) :: nx
      real(kind=8), dimension(nx), intent(out) :: x_q, w_q
      real(kind=8), dimension(6,24) :: x, a
      integer :: j, n
  
      ! Initialize x and a arrays (you need to fill these with the correct values)
      ! ...
  
      n = nx / 4
      nx = n * 4
      if (n < 1 .or. n > 6) then
          write(*, '(A,I0)') "N/4 out of range : ", n
          return
      end if
  
      do j = 1, nx
          x_q(j) = x(n, j)
          w_q(j) = a(n, j)
      end do
  end subroutine laguerreq
  
  