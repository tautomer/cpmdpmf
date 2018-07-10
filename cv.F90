function proj(r)
    implicit none

    real*8, intent(in) :: r(3, 3)
    real*8 :: xy(3), xz(3), lxz, dotxyxz
    real*8 :: proj

    lxz = 0
    dotxyxz = 0
    xy = r(1, :) - r(2, :)
    xz = r(2, :) - r(3, :)
    lxz = lxz + xz(1) ** 2
    lxz = lxz + xz(2) ** 2
    lxz = lxz + xz(3) ** 2
    lxz = dsqrt(lxz)
    dotxyxz = dotxyxz + xy(1) * xz(1)
    dotxyxz = dotxyxz + xy(2) * xz(2)
    dotxyxz = dotxyxz + xy(3) * xz(3)
    proj = dotxyxz / lxz

    return
end function

function diff(r)
    implicit none

    real*8, intent(in) :: r(3, 3)
    real*8 :: xy(3), yz(3), lxy, lyz
    real*8 :: diff

    lxy = 0
    lyz = 0
    xy = r(1, :) - r(2, :)
    yz = r(2, :) - r(3, :)
    lxy = lxy + xy(1) ** 2
    lxy = lxy + xy(2) ** 2
    lxy = lxy + xy(3) ** 2
    lyz = lyz + yz(1) ** 2
    lyz = lyz + yz(2) ** 2
    lyz = lyz + yz(3) ** 2
    lxy = dsqrt(lxy)
    lyz = dsqrt(lyz)
    diff = lxy - lyz

    return
end function
