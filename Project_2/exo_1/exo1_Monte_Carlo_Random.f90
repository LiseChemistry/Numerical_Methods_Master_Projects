program monte_carlo
implicit none
real(8) :: r
real(8) :: sum
real(8), parameter :: l = 0.0d0 ! orbital s
real(8), parameter :: n = 3.0d0 ! ns = 3n
real(8), parameter :: a_o = 1.0d0 ! atomic units
real(8), parameter :: pi = 3.141592653589793d0
real(8), parameter :: Z = 1.0d0 ! Charge
integer, parameter :: M = 1000000 ! Number of sampless
integer :: i
real(8) :: V = 61.0d0
real(8) :: phi_3s
real(8) :: moy_r
real(8) :: Random, moy_f
 
sum = 0.0d0
phi_3s = 0.0d0
moy_f = 0.0d0
moy_r = 0.0d0

! To generate the random nubers
do i = 1, M
    Random = rand() * 61.0d0
    r = Random
    phi_3s =(1.0d0/ sqrt (4.0d0*pi)) * (2.0d0 / (81.0d0 * sqrt (3.0d0))) * &
    (27.0d0 - (18.0d0 * r) + 2.0d0*(r**2.0d0)) * &
    exp(-(r/3.0d0))
    moy_r = 4.0d0 * pi * (((phi_3s)**2.0d0) * (r**3.0d0))
    print *, moy_r , i
    sum = sum + moy_r
end do
    ! To calculate the weighted mean value
    moy_f = (1.0d0/ M) * sum * 61.0d0

print *, "The mean distance ⟨r⟩ 3s is :", moy_f

end program monte_carlo
