subroutine transpon(transpon_tensor, tensor)
    !input: tensor
    !output: tensor^T
    real*8 :: transpon_tensor(3,3), tensor(3,3)
    integer :: i,j
    do i=1,3
        do j=1,3
            transpon_tensor(j,i)=tensor(i,j)
        end do
    end do
end subroutine transpon    