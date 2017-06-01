!$openad XXX Template ad_template.f
subroutine head(s,zres) 
	real(kind=8), dimension(3) :: X
	real(kind=8), dimension(3) :: Xnp1_res
	real(kind=8), intent(in) :: s
	real(kind=8), intent(out) :: zres
!	use Lorenz63
	X(1) = -2.4d0
	X(2) = -3.7d0
	X(3) = 14.98d0
	do t = 1, 1000, 1
		call Xnp1(X,Xnp1_res,s)
	end do


!$openad INDEPENDENT(s) 
call z(Xnp1_res,zres,s) 
!$openad DEPENDENT(zres)
end subroutine

!$openad XXX Template ad_template.f
subroutine z(X,z_res,r)

    implicit none
	real(kind=8), dimension(3) :: X
    real(kind=8), dimension(3):: k1, k2, k3, k4
	real(kind=8), dimension(3) :: Xnp1_res
    real(kind=8), dimension(3):: ddt
	integer :: i
	real(kind=8) :: dt
	real(kind=8), intent(in) :: r
	real(kind=8), intent(out):: z_res	
	dt = 0.005d0
		
    call dXdt(X,ddt,r)
    k1 = dt*ddt
    call dXdt(X+0.5d0*k1,ddt,r)
    k2 = dt*ddt
    call dXdt(X+0.5d0*k2,ddt,r)
    k3 = dt*ddt
    call dXdt(X+k3,ddt,r)
    k4 = dt*ddt
    
    Xnp1_res = X + 1.d0/6.d0*k1 + &
               1.d0/3.d0*k2 + 1.d0/3.d0*k3 + &
                1.d0/6.d0*k4   


	z_res = Xnp1_res(3)

end subroutine z

!$openad XXX Template ad_template.f
subroutine Xnp1(X,Xnp1_res,r)

        
    implicit none
	real(kind=8), intent(in), dimension(3) :: X
    real(kind=8), dimension(3):: k1, k2, k3, k4
	real(kind=8), intent(out), dimension(3) :: Xnp1_res
    real(kind=8), dimension(3):: ddt
	integer :: i
	real(kind=8) :: dt
	real(kind=8), intent(in) :: r
	
	dt = 0.005d0
		
    call dXdt(X,ddt,r)
    k1 = dt*ddt
    call dXdt(X+0.5d0*k1,ddt,r)
    k2 = dt*ddt
    call dXdt(X+0.5d0*k2,ddt,r)
    k3 = dt*ddt
    call dXdt(X+k3,ddt,r)
    k4 = dt*ddt
    
    Xnp1_res = X + 1.d0/6.d0*k1 + &
               1.d0/3.d0*k2 + 1.d0/3.d0*k3 + &
                1.d0/6.d0*k4   



end subroutine Xnp1
subroutine sys_params(sigma, b)
	
	implicit none
	real(kind=8), intent(out) :: sigma, b
	sigma = 10.d0
	b = 8.d0/3.d0

end subroutine sys_params

!$openad XXX Template ad_template.f
subroutine dXdt(X,dXdt_res,r)
    implicit none
	real(kind=8), intent(in), dimension(3) :: X
	real(kind=8), intent(in) :: r
	real(kind=8), intent(out), dimension(3) :: dXdt_res
	real(kind=8) :: sigma, b
	integer :: i
	real(kind=8) :: dt
    
    dt = 0.005d0
	call sys_params(sigma,b) 	
	dXdt_res(1) = -sigma*X(1) + sigma*X(2)
	dXdt_res(2) = X(1)*(r-X(3)) - X(2)
	dXdt_res(3) = X(1)*X(2) - b*X(3)  
        
end subroutine dXdt

