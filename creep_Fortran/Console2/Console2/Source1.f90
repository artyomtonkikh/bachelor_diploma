program prog
    real*8 :: x(10)
    integer :: i, k
    x(1)=1
    do i=2,10
       x(i)=x(i-1)+i 
    end do
    write (*,*) x
    read (*,*) k
    x=k*x
    write (*,*) x
    read (*,*) k
end