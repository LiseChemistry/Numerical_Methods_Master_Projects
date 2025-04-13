module AnalyticFunctionModule
    implicit none
    real(kind=8), parameter:: pi = 3.1415926535897932d0
    type, public :: AnalyticFunction
        private
        real(kind=8), dimension(:), allocatable :: x_grid
        real(kind=8), dimension(:), allocatable :: f_grid
        integer(kind=8) :: nb_point
        real(kind=8) :: upper_bound
    contains
        procedure, public :: compute, integrate, initialize, get_nb_point, get_upper_bound
    end type AnalyticFunction
contains
    subroutine initialize(this, x_grid, f_grid, nb_point, upper_bound)
        class(AnalyticFunction), intent(inout) :: this
        real(kind=8), dimension(:), allocatable, intent(in) :: x_grid, f_grid
        integer(kind=8), intent(in) :: nb_point
        real(kind=8), intent(in) :: upper_bound
        this%x_grid = x_grid
        this%f_grid = f_grid
        this%nb_point = nb_point
        this%upper_bound = upper_bound
    end subroutine

    function get_nb_point(this)
        class(AnalyticFunction), intent(in) :: this
        integer(kind=8) :: get_nb_point
        get_nb_point = this%nb_point
    end function get_nb_point

    function get_upper_bound(this)
        class(AnalyticFunction), intent(in) :: this
        real(kind=8) :: get_upper_bound
        get_upper_bound = this%upper_bound
    end function get_upper_bound


    function compute(this, x)
        class(AnalyticFunction), intent(in) :: this
        real(kind=8), intent(in) :: x
        integer(kind=8) :: i
        real(kind=8) :: compute
        real(kind=8) :: h

        compute = 0.0d0
        h = (this%upper_bound) / (this%nb_point-1)

        do i = 1, this%nb_point
            compute = compute + (-1)**(i-1) * this%f_grid(i) * sin(x * pi / h) /(x- this%x_grid(i))
        end do
        compute = compute * h / pi

    end function

    function integrate(this)
        class(AnalyticFunction), intent(in) :: this
        integer(kind=8) :: i
        real(kind=8) :: integrate
        real(kind=8) :: h

        integrate = 0.0d0
        h = this%upper_bound / (this%nb_point-1)

        do i = 1, this%nb_point
            integrate = integrate + (-1)**(i-1) * this%f_grid(i) * pi * cos(pi * this%x_grid(i) / h)
        end do
        integrate = integrate * h / pi

    end function

end module AnalyticFunctionModule

module F_IntegrandModule
    implicit none
    type, public :: F_Integrand
        private
        real(kind=8) :: omega
        real(kind=8) :: alpha
        real(kind=8) :: T
    contains
        procedure, public :: initialize, compute
    end type F_Integrand

contains
    subroutine initialize(this, omega, alpha, T)
        class(F_Integrand), intent(inout) :: this
        real(kind=8), intent(in) :: omega
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: T
        this%omega = omega
        this%alpha = alpha
        this%T = T
    end subroutine

    function compute(this, t)
        class(F_Integrand), intent(in) :: this
        real(kind=8), intent(in) :: t
        real(kind=8) :: compute
        compute = sin(this%omega * t) * exp(-this%alpha * ((t-this%T)**2))
        compute = compute * compute
    end function

end module F_IntegrandModule

!#############################################################################################################################################
!GridPrinterModule
!
! procedure :
! print_grid : print vector x, f in a file
!##############################################################################################################################################

module GridPrinterModule
    implicit none
contains
  subroutine print_grid(file_name, x_grid, y_grid, grid_size)

        character(len=:), allocatable, intent(in) :: file_name
        integer(kind=8), intent(in) :: grid_size
        real(kind=8), dimension (grid_size), intent(in) :: x_grid
        real(kind=8), dimension (grid_size), intent(in) :: y_grid

        integer(kind=8) :: i, unit_number
        unit_number = 1

        open(unit = unit_number, file=file_name, form='formatted')

        write(unit_number, *) '#' , ' x', ' y'
        do i =1, grid_size
            write(unit_number, *) x_grid(i) , y_grid(i)
        end do
        close(unit = unit_number)
    end subroutine

end module GridPrinterModule

module SimpsonIntegrationModule
    use F_IntegrandModule, only : F_Integrand
    implicit none
    type, public :: SimpsonIntegration
        private
        integer(kind=8) :: nb_point
    contains
        procedure, public :: initialize, integrate
    end type SimpsonIntegration

contains
    subroutine initialize(this, nb_point)
        class(SimpsonIntegration), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_point
        integer(kind=8), parameter:: two = 2

        if(nb_point.le.two) then
            print* , 'Simpson Integration needs at least 3 points'
        end if
        if(modulo(nb_point, two) == 0) then
            print* , 'This Simpson Integration implementation is designed for an odd number of point'
        end if

        this%nb_point = nb_point
    end subroutine

    function integrate(this, integrand, lower_bound, upper_bound)
        class(SimpsonIntegration), intent(in) :: this
        class(F_Integrand), intent(in) :: integrand
        real(kind=8), intent(in) :: lower_bound, upper_bound
        real(kind=8) :: integrate
        real(kind=8) :: step
        integer(kind=8) :: i

        step = (upper_bound - lower_bound) / real(this%nb_point - 1, 8)

        integrate = 0.5d0 * integrand%compute(lower_bound)

        do i = 2, this%nb_point - 3, 2
            integrate = integrate + integrand%compute(lower_bound + i *step)
        end do

        do i = 1, this%nb_point - 2, 2
            integrate = integrate + 2.0d0 * integrand%compute(lower_bound + i *step)
        end do

        integrate = integrate + 0.5d0 * integrand%compute(upper_bound)
        integrate = integrate * step * 2.0d0 / 3.0d0

    end function


end module SimpsonIntegrationModule

module TrapezoidalIntegrationModule
    use F_IntegrandModule, only : F_Integrand
    implicit none
    type, public :: TrapezoidalIntegration
        private
        integer(kind=8) :: nb_point
    contains
        procedure, public :: initialize, integrate
    end type TrapezoidalIntegration

contains
    subroutine initialize(this, nb_point)
        class(TrapezoidalIntegration), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_point
        if(nb_point.le.1) then
            print* , 'Trapezoidal Integration needs at least 2 points'
        end if
        this%nb_point = nb_point
    end subroutine

    function integrate(this, integrand, lower_bound, upper_bound)
        class(TrapezoidalIntegration), intent(in) :: this
        class(F_Integrand), intent(in) :: integrand
        real(kind=8), intent(in) :: lower_bound, upper_bound
        real(kind=8) :: integrate
        real(kind=8) :: step
        integer(kind=8) :: i

        step = (upper_bound - lower_bound) / real(this%nb_point - 1, 8)

        integrate = 0.5d0 * integrand%compute(lower_bound)

        do i = 1, this%nb_point - 2
            integrate = integrate + integrand%compute(lower_bound + i *step)
        end do

        integrate = integrate + 0.5d0 * integrand%compute(upper_bound)
        integrate = integrate * step

    end function

end module TrapezoidalIntegrationModule

program exercice_3
    use F_IntegrandModule, only : F_Integrand
    use TrapezoidalIntegrationModule
    use SimpsonIntegrationModule
    use AnalyticFunctionModule, only : AnalyticFunction, pi
    use GridPrinterModule
    implicit none

    type(F_Integrand) :: integrand
    real(kind=8), parameter :: omega = 2.0d0
    real(kind=8), parameter :: alpha = 0.05
    real(kind=8), parameter :: T =  5.0d0 * pi
    real(kind=8), parameter :: upper_bound =  2.0d0*T

    real(kind=8), parameter :: eps = 1.0d-16
    real(kind=8) :: integral_n,integral_n_plus_1, iteration_difference
    real(kind=8), parameter :: threshold = 1.0d-8
    integer(kind=8) :: nb_point
    integer(kind=8) i, j
    real(kind=8) :: x_space_step
    integer(kind=8), parameter :: grid_size_print = 10000
    real(kind=8), dimension(grid_size_print) :: x_grid_print, f_analytic_grid_print, f_integrand_grid_print

    integer(kind=8), dimension(2) :: grid_size_analytic
    real(kind=8), dimension(:), allocatable :: x_grid, f_grid
    real(kind=8) :: h
    type(AnalyticFunction) :: analytic_function

    character(len=:), allocatable :: print_file
    type(TrapezoidalIntegration) trapezoidal_integration
    type(SimpsonIntegration) simpson_integration
    character(len=10) :: chaine_grid_size


    !Question 1
    call integrand%initialize(omega, alpha, T)

    nb_point = 1
    iteration_difference = 1.0d0
    do while((iteration_difference > threshold).or.(abs(integral_n)<eps))

        nb_point = nb_point + 1

        call trapezoidal_integration%initialize(nb_point)
        integral_n = trapezoidal_integration%integrate(integrand, 0.0d0, upper_bound)

        call trapezoidal_integration%initialize(nb_point +1)
        integral_n_plus_1 = trapezoidal_integration%integrate(integrand, 0.0d0, upper_bound)

        iteration_difference = abs(integral_n_plus_1 - integral_n)

    end do

    print*, 'Results trapeze'
    print*, 'Number of points n : ', nb_point
    print*, 'I(n+1) - I(n) : ' , iteration_difference
    print*, 'I(n) : ', integral_n


    nb_point = 1
    iteration_difference = 1.0d0
    do while((iteration_difference > threshold).or.(abs(integral_n)<eps))

        nb_point = nb_point + 2

        call simpson_integration%initialize(nb_point)
        integral_n = simpson_integration%integrate(integrand, 0.0d0, upper_bound)

        call simpson_integration%initialize(nb_point + 2)
        integral_n_plus_1 = simpson_integration%integrate(integrand, 0.0d0, upper_bound)

        iteration_difference = abs(integral_n_plus_1 - integral_n)

    end do

    print*, 'Results Simpson'
    print*, 'Number of points n : ', nb_point
    print*, 'I(n+2) - I(n) : ' , iteration_difference
    print*, 'I(n) : ', integral_n


    !Question 2
    grid_size_analytic(1) = 20 ! test avec 20 points
    grid_size_analytic(2) = 60 ! test avec 60 points

    do j = 1, 2

        h = upper_bound /real(grid_size_analytic(j)-1,8)

        allocate(x_grid(grid_size_analytic(j)))
        allocate(f_grid(grid_size_analytic(j)))

        do i =1, grid_size_analytic(j)
            x_grid(i) = real((i-1),8) * h
            f_grid(i) = integrand%compute(x_grid(i))
        end do

        call analytic_function%initialize(x_grid, f_grid, grid_size_analytic(j), upper_bound)

        x_space_step =  upper_bound / real(grid_size_print - 1 , 8)

        print*, grid_size_analytic(j), 'points' , ', Results : ',  analytic_function%integrate()

        do i =1, grid_size_print
            x_grid_print(i) = x_space_step * real((i -1),8)
            f_analytic_grid_print(i) = analytic_function%compute(x_grid_print(i))
            f_integrand_grid_print(i) = integrand%compute(x_grid_print(i))
        end do


        write(chaine_grid_size, '(I0)') grid_size_analytic(j)

        print_file  = "..\Results\analytic_function_" // trim(chaine_grid_size) // ".dat"
        call print_grid(print_file , x_grid_print, f_analytic_grid_print, grid_size_print)

        print_file  = "..\Results\integrand_" // trim(chaine_grid_size) // ".dat"
        call print_grid(print_file , x_grid_print, f_integrand_grid_print, grid_size_print)

        deallocate(x_grid)
        deallocate(f_grid)

    end do

end program
