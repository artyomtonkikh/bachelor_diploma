function trace(tensor)
    !input: tensor
    !output: tr(tensor)
    real*8 :: tensor(3,3)
    real*8 :: trace
    integer :: i
    trace=0.0d0
    do i=1,3
        trace=trace+tensor(i,i)
    end do
end function trace    