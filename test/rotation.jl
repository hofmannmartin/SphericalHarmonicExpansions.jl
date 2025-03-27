############################################################################
## Testing the rotation and reflection of spherical harmonic coefficients ##
############################################################################
using Rotations # rotated positions

@testset "Rotation and reflection" begin

  # test different rotations
  @testset "Rotation" begin
    ## test if expansion with rotated coefficients yields same values  
    ## as the expansion with rotated positions
    c = SphericalHarmonicCoefficients(randn(16))
    p = fastfunc(c)

    # random rotation of the coefficients
    ang = randn(3) # random Euler angles (ZYZ convention)
    c_rot = rotation(c, ang)
    p_rot = fastfunc(c_rot)

    # evaluate functions
    pos = randn(3) # position in the rotated coordinate system
    pos_oldCS = transpose(RotZYZ(ang...)) * pos # position in the old coordinate system
    val = p(pos_oldCS)
    val_rot = p_rot(pos)

    # test
    @test isapprox(val, val_rot)
  end

  # test the reflection at the origin and different planes
  @testset "Point reflection" begin
    val = randn(16)
    c = SphericalHarmonicCoefficients(val)

    ## 1. Reflection at the origin
    c_ref = pointReflection(c)
    # correct solution
    val_solution = [val[1], -val[2:4]..., val[5:9]..., -val[10:16]...]
    c_solution = SphericalHarmonicCoefficients(val_solution)
    # test
    @test isapprox(c_ref, c_solution)

    ## 2. Second reflection at the origin should give the original coefficients
    c_ref = pointReflection(c_ref)
    @test isapprox(c, c_ref)
  end

  # Reflection at different planes
  @testset "Plane reflection" begin
    val = randn(9)
    c = SphericalHarmonicCoefficients(val)

    # 1. Reflection at the xz-plane (y -> -y, i.e., apply a point reflection and a subsequent rotation around the y-axis by pi)
    c_rot = rotation(pointReflection(c), [0,pi,0]);
    # correct solution
    val_solution = [val[1], -val[2], val[3:4]..., -val[5:6]..., val[7:9]...]
    c_solution = SphericalHarmonicCoefficients(val_solution)
    # test
    @test isapprox(c_rot, c_solution)

    # 2. Reflection at the xy-plane (z -> -z, i.e., apply a point reflection and a subsequent rotation around the z-axis by pi)
    c_rot = rotation(pointReflection(c), [pi,0,0]);
    # correct solution
    val_solution = [val[1:2]..., -val[3], val[4:5]..., -val[6], val[7], -val[8], val[9]]
    c_solution = SphericalHarmonicCoefficients(val_solution)
    # test
    @test isapprox(c_rot, c_solution)

    # 3. Reflection at the yz-plane (x -> -x, i.e., apply a point reflection and a subsequent rotation around the x-axis by pi)
    c_rot = rotation(pointReflection(c), [pi,pi,0]);
    # correct solution
    val_solution = [val[1], val[2:3]..., -val[4:5]..., val[6:7]..., -val[8], val[9]]
    c_solution = SphericalHarmonicCoefficients(val_solution)
    # test
    @test isapprox(c_rot, c_solution)
  end
end