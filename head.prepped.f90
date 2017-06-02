subroutine head(s,zres)
	use Lorenz63
	implicit none 
	real(kind=8), dimension(3) :: X
	real(kind=8), dimension(3) :: Xnp1_res
	real(kind=8), intent(in) :: s
	real(kind=8), intent(out) :: zres
	integer :: t



!$openad INDEPENDENT(s) 
	X(1) = -2.4d0
	X(2) = -3.7d0
	X(3) = 14.98d0
	do t = 1, 1000, 1
		call Xnp1(X,Xnp1_res,s)
		X = Xnp1_res
	end do

	zres = Xnp1_res(3) 

!$openad DEPENDENT(zres)
end subroutine


