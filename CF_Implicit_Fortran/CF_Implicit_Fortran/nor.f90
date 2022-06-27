real*8 function nor(tensor)
    !input: tensor
    !output: ||tensor||
    real*8 :: tensor(3,3)
    real*8 :: norma
    integer :: i,j
    nor=0.0d0
    do i=1,3
        do j=1,3
          nor=nor+tensor(i,j)**2.0d0 
        end do
    end do
    nor=sqrt(nor)
end function nor 