Module global_var

    use size 

    IMPLICIT NONE

    DOUBLE PRECISION :: h(0:nx+1)
    DOUBLE PRECISION :: A(0:nx+1)
    DOUBLE PRECISION :: hn1(0:nx+1)
    DOUBLE PRECISION :: An1(0:nx+1)
    DOUBLE PRECISION :: hn2(0:nx+1)
    DOUBLE PRECISION :: An2(0:nx+1)
    DOUBLE PRECISION :: Pp(0:nx+1)
    DOUBLE PRECISION :: P(0:nx+1)
    DOUBLE PRECISION :: sig11(0:nx+1)

    ! Muphi
    DOUBLE PRECISION :: shear_I(0:nx+1)
    DOUBLE PRECISION :: mu_I(0:nx+1)
    DOUBLE PRECISION :: Inertial(0:nx+1)

end Module global_var