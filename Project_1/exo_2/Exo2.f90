!######################################################################################################
! type ModelPotentialCalculator
! data members :
! alpha : parameter which describe model potential, see formula in the report
! procedure members :
! initialize : it is called to set data member before doing a computation
! compute : return the potential value at phi
!######################################################################################################

module ModelPotentialCalculatorModule
    implicit none
    type, public :: ModelPotentialCalculator
        private
        real(kind=8) :: alpha
    contains
        procedure, public :: initialize, compute
    end type ModelPotentialCalculator

contains
    subroutine initialize(this, alpha)
        class(ModelPotentialCalculator), intent(inout) :: this
        real(kind=8), intent(in) :: alpha
        this%alpha = alpha
    end subroutine

    function compute(this, phi)
        class(ModelPotentialCalculator), intent(in) :: this
        real(kind=8), intent(in) :: phi
        real(kind=8) compute
        compute = exp(- this%alpha * cos(6.0d0 * phi))
    end function
end module ModelPotentialCalculatorModule

!#############################################################################################################################################
!GridPrinterModule
!
! procedure :
! print_grid : print vector x, f in a file
! print_grid_with_error : print vector x, f, error in a file
!##############################################################################################################################################

module GridPrinterModule
    implicit none
contains
    subroutine print_grid(file_name, x_grid, f_grid, grid_size)

        character(len=:), allocatable, intent(in) :: file_name
        integer, intent(in) :: grid_size
        real(kind=8), dimension (grid_size), intent(in) :: x_grid
        real(kind=8), dimension (grid_size), intent(in) :: f_grid
        integer :: i, unit_number
        unit_number = 1

        open(unit = unit_number, file=file_name, form='formatted')

        write(unit_number, *) '#' , ' x', ' f'
        do i =1, grid_size
            write(unit_number, *) x_grid(i) , f_grid(i)
        end do
    end subroutine

    subroutine print_grid_with_error(file_name, x_grid, f_grid, error_grid, grid_size)

        character(len=:), allocatable, intent(in) :: file_name
        integer, intent(in) :: grid_size
        real(kind=8), dimension (grid_size), intent(in) :: x_grid
        real(kind=8), dimension (grid_size), intent(in) :: f_grid
        real(kind=8), dimension (grid_size), intent(in) :: error_grid
        integer :: i, unit_number
        unit_number = 1

        open(unit = unit_number, file=file_name, form='formatted')

        write(unit_number, *) '#' , ' x', ' f', ' error'
        do i =1, grid_size
            write(unit_number, *) x_grid(i) , f_grid(i), error_grid(i)
        end do
    end subroutine

end module GridPrinterModule


!####################################################################################
! type C2H5PotentialTrigonometricInterpolator
! data members :
! grid_size : number of points used
! coefficient : coefficient of the cosinus basis functions
!
! procedure members :
! build_phi_grid : it is used to build phi grid according to Gauss-Tch�bytcheff quadrature
! initialize : it is called to set data member before doing an interpolation
!####################################################################################
!   ===+=========+=========+=========+=========+=========+=========+==
    Subroutine Minv(A,nA,dimA,L,Wk,tol)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Matrix inversion and solution of linear systems
!     A(dimA,1) : matrix to be inverted
!     dimA : dimension declared in the calling sequence
!     nA : dimension of the system
!     L    : work array of dimension 2*n
!     Wk   :  "     "   "      "     " .On return, Wk(1) contains the
!            determinant.
!     Tol : tolerence on pivoting (to be set to 0d0 in the call)
!     ===+=========+=========+=========+=========+=========+=========+==
      Implicit Real*8(a-h,o-z)
      Integer dimA
      Dimension A(*),L(*),Wk(*)
!        search for largest element
!
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
!
!        interchange rows
!
      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-dima
      do 30 i=1,nA
      ki=ki+dima
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold
!
!        interchange columns
!
   35 i=l(nA+k)
      if(i-k) 45,45,38
   38 jp=dima*(i-1)
      do 40 j=1,nA
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold
!
!        divide column by
!        contained in biga)
!
   45 if(abs(biga).le.tol) then
          wk(1)=0.
          return
      endif
      do 55 i=1,nA
      if(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue
!
!        reduce matrix
!
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
!
!        divide row by pivot
!
      kj=k-dima
      do 75 j=1,nA
      kj=kj+dima
      if(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue
!
!        product of pivots
      d=d*biga
!
!        replace pivot by reciprocal
      a(kk)=1.d0/biga
   80 continue
!
!        final row and column interchange
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
!
  150 Wk(1)=d
      Return
      
end Subroutine

module C2H5PotentialTrigonometricInterpolatorModule
    implicit none
    real(kind=8), parameter:: pi = 3.1415926535897932d0
    type, public :: C2H5PotentialTrigonometricInterpolator
        private
        integer :: grid_size
        real(kind=8), dimension(:), allocatable :: coefficient
    contains
      procedure, public ::  initialize, interpolate
    end type C2H5PotentialTrigonometricInterpolator

contains

     function build_phi_grid(grid_size)
        integer :: grid_size
        real(kind=8), dimension(grid_size) :: build_phi_grid
        integer :: i

        do i=1, grid_size
            build_phi_grid(i) = (2.0d0 * real(i-1, kind=8) + 1.0d0) * pi / (12.0d0 * real(grid_size, kind=8))
        end do
    end function

    subroutine initialize(this, grid_size, phi_grid, v_grid)
        class(C2H5PotentialTrigonometricInterpolator), intent(inout) :: this
        integer, intent(in) :: grid_size
        real(kind=8), dimension(grid_size), intent(in) :: phi_grid, v_grid
        real(kind=8), dimension(grid_size, grid_size) :: r_matrix
        real(kind=8), dimension(grid_size, grid_size) :: inverse_r_matrix
        real(kind=8), dimension(2 * grid_size) :: Wk
        integer, dimension(2 * grid_size) :: L
        real(kind=8), dimension(grid_size) :: inverse_r_times_v
        real(kind=8), parameter :: zero_threshold = 1.0d-16
        integer :: i, j

        this%grid_size = grid_size

        do i=1, grid_size
            do j=1, grid_size
                r_matrix(i,j)= cos(6.0d0 * real((j-1), kind=8) * phi_grid(i))
            end do
        end do
        inverse_r_matrix = r_matrix

        call Minv(inverse_r_matrix, grid_size, grid_size, L, Wk, 0.0d0)

        !Check if the matrix is invertible (Wk(1) contains determinant of the matrix)
        if (abs(Wk(1)).le.zero_threshold) then
            print *, 'Initialization failed as R matrix is not invertible'
            return
        end if

        inverse_r_times_v = matmul(inverse_r_matrix, v_grid)
        this%coefficient = inverse_r_times_v
    end subroutine

    function interpolate(this, phi)
        class(C2H5PotentialTrigonometricInterpolator), intent(in) :: this
        real(kind=8), intent(in) :: phi
        real(kind=8) :: interpolate
        integer :: i
        interpolate = 0.0d0

        do i = 1, this%grid_size
            interpolate = interpolate + this%coefficient(i) * cos(6.0d0 * real((i-1), kind=8) * phi)
        end do
    end function

end module

program project_1_2
    use C2H5PotentialTrigonometricInterpolatorModule
    use ModelPotentialCalculatorModule
    use GridPrinterModule
    implicit none

    type(C2H5PotentialTrigonometricInterpolator) :: trigonometric_interpolator
    type(ModelPotentialCalculator) :: model_potential_calculator
    real(kind=8), dimension(:), allocatable :: phi_grid, v_grid
    real(kind=8), parameter :: alpha = 0.2d0
    real(kind=8), parameter :: threshold = 1.0d-5
    real(kind=8) :: stddev, intermediate_phi, interpolation_error, x_space_step
    integer :: i, grid_size, intermediate_grid_size
    real(kind=8), parameter :: x_grid_minimum = 0.0d0, x_grid_maximum = 2.0d0 * pi
    integer, parameter :: grid_size_output = 500

    real(kind=8), dimension(grid_size_output) :: x_grid, potential_f_grid, trigonometric_f_grid, &
                                                 error_grid

    character(len=:), allocatable :: model_potential_file, trigonometric_file
    character(len=2) :: chaine_grid_size

    grid_size = 1
    stddev = 1.0d0
    x_space_step = (x_grid_maximum - x_grid_minimum) / real(grid_size_output - 1 , 8)

    call model_potential_calculator%initialize(alpha)

    !################################################################################
    !
    ! Find the minimum number of points that satisfies the criteria of the exercise
    !
    !################################################################################

    do while(stddev > threshold)

        grid_size = grid_size + 1
        phi_grid = build_phi_grid(grid_size)
        allocate(v_grid(grid_size))
        do i = 1, grid_size
            v_grid(i) = model_potential_calculator%compute(phi_grid(i))
        end do

        call trigonometric_interpolator%initialize(grid_size, phi_grid, v_grid)

        intermediate_grid_size = grid_size - 1
    
        print *, 'Trigonometric error check at grid point for grid size =', grid_size
        do i =1, grid_size
            print*,  abs(trigonometric_interpolator%interpolate(phi_grid(i))- &
            model_potential_calculator%compute(phi_grid(i)))
        end do

        stddev = 0.0d0
        do i =1, intermediate_grid_size
            intermediate_phi = (phi_grid(i+1) + phi_grid(i)) / 2.0d0
            interpolation_error = trigonometric_interpolator%interpolate(intermediate_phi) - &
            model_potential_calculator%compute(intermediate_phi)
            stddev = stddev + interpolation_error **2
        end do

        stddev = stddev / intermediate_grid_size
        stddev = sqrt(stddev)
        deallocate(v_grid)

        do i =1, grid_size_output
            x_grid(i) = x_grid_minimum + x_space_step * real((i -1),8)
            potential_f_grid(i) =  model_potential_calculator%compute(x_grid(i))
            trigonometric_f_grid(i) = trigonometric_interpolator%interpolate(x_grid(i))

            error_grid(i) = abs(potential_f_grid(i) - trigonometric_f_grid(i))
        end do

        write(chaine_grid_size, '(I0)') grid_size
        trigonometric_file = "..\\Results\trigonometric" // trim(chaine_grid_size) // ".dat"
        call print_grid_with_error(trigonometric_file, x_grid, trigonometric_f_grid, error_grid, grid_size_output)
        print *, 'Number of point : ', grid_size
        print *, 'Standard deviation : ', stddev
    end do

    !#####################################################
    !
    ! Print result in files
    !
    !#####################################################
    model_potential_file = "..\\Results\model_potential_2.dat"

    do i =1, grid_size_output
        x_grid(i) = x_grid_minimum + x_space_step * real((i -1),8)
        potential_f_grid(i) =  model_potential_calculator%compute(x_grid(i))
    end do
    call print_grid(model_potential_file, x_grid, potential_f_grid, grid_size_output)
    call system ('gnuplot -p potentiel.plt')
end program


