program driver
  use OAD_active
  use OAD_rev
  implicit none 
  external head
  type(active) :: x, y
  x%v=2.8D1
  y%d=1.0D0
  our_rev_mode%tape=.TRUE.
  call head(x,y)
  print *, 'driver running for x =',x%v
  print *, '            yields y =',y%v,' dy/dx =',x%d
  print *, '    dz^1000/dz = ',1.D0*cos(x%v)-x%d
end program driver
