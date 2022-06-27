subroutine dev(dev_tensor, tensor)
    !input: tensor
    !output: dev_tensor
    real*8 :: tr, trace
    real*8 :: dev_tensor(3,3), tensor(3,3)
    real*8 :: eye(3,3) !единичный тензор
    integer :: i
    eye=0.0d0
    do i=1,3
        eye(i,i)=1.0d0
    end do
    tr=trace(tensor)
    dev_tensor=tensor-tr*1.0d0/3.0d0*eye
end subroutine dev    