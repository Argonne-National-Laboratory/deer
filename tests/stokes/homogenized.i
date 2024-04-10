[Mesh]
  second_order = true
  [mesh]
    type = CartesianMeshGenerator
    dim = 3
    dx = '1 1'
    dy = '1 1'
    dz = '1 1'
    ix = '8 8'
    iy = '8 8'
    iz = '8 8'
    subdomain_id = '1 2 3 4 5 6 7 8'
  []
[]

[Variables]
  [p]
    order = FIRST
    family = LAGRANGE
  []
  [u]
    order = SECOND
    family = LAGRANGE_VEC
  []
[]

[AuxVariables]
    [deviatoric_stress_xx]
        family = MONOMIAL
        order = CONSTANT
    []
    [deviatoric_stress_yy]
        family = MONOMIAL
        order = CONSTANT
    []
    [deviatoric_stress_zz]
        family = MONOMIAL
        order = CONSTANT
    []
    [deviatoric_stress_xy]
        family = MONOMIAL
        order = CONSTANT
    []
    [deviatoric_stress_yz]
        family = MONOMIAL
        order = CONSTANT
    []
    [deviatoric_stress_xz]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_xx]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_yy]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_zz]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_xy]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_yz]
        family = MONOMIAL
        order = CONSTANT
    []
    [stress_xz]
        family = MONOMIAL
        order = CONSTANT
    []
[]

[AuxKernels]
    [deviatoric_stress_xx]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_xx
        property = 'deviatoric_stress'
        i = 0
        j = 0
    []
    [deviatoric_stress_yy]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_yy
        property = 'deviatoric_stress'
        i = 1
        j = 1
    []
    [deviatoric_stress_zz]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_zz
        property = 'deviatoric_stress'
        i = 2
        j = 2
    []
    [deviatoric_stress_xy]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_xy
        property = 'deviatoric_stress'
        i = 0
        j = 1
    []
    [deviatoric_stress_yz]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_yz
        property = 'deviatoric_stress'
        i = 1
        j = 2
    []
    [deviatoric_stress_xz]
        type = ADMaterialRankTwoTensorAux
        variable = deviatoric_stress_xz
        property = 'deviatoric_stress'
        i = 0
        j = 2
    []
    [stress_xx]
      type = ParsedAux
      variable = stress_xx
      coupled_variables = 'p deviatoric_stress_xx'
      expression = 'deviatoric_stress_xx - p'
    []
    [stress_yy]
      type = ParsedAux
      variable = stress_yy
      coupled_variables = 'p deviatoric_stress_yy'
      expression = 'deviatoric_stress_yy - p'
    []
    [stress_zz]
      type = ParsedAux
      variable = stress_zz
      coupled_variables = 'p deviatoric_stress_zz'
      expression = 'deviatoric_stress_zz - p'
    []
    [stress_xy]
      type = ParsedAux
      variable = stress_xy
      coupled_variables = 'deviatoric_stress_xy'
      expression = 'deviatoric_stress_xy'
    []
    [stress_xz]
      type = ParsedAux
      variable = stress_xz
      coupled_variables = 'deviatoric_stress_xz'
      expression = 'deviatoric_stress_xz'
    []
    [stress_yz]
      type = ParsedAux
      variable = stress_yz
      coupled_variables = 'deviatoric_stress_yz'
      expression = 'deviatoric_stress_yz'
    []
[]

[ICs]
  [u]
    type = VectorFunctionIC
    variable = u
    function_x = '1.0 * x'
    function_y = '1.0 * y'
    function_z = '1.0 * z'
  []
[]

[Kernels]
  [equil]
    type = ADStokesStressDivergence
    variable = u
  []
  [pressure]
    type = ADStokesPressure
    variable = u
    pressure = p
  []
  [incompressible]
    type = ADStokesIncompressibility
    variable = p
    velocity = u
  []
[]

[Materials]
  [strain]
    type = StokesStrainRate
    velocity = u
  []
  [stress1]
    type = StokesLinearViscous
    mu = 1
    block = '1 4 6 7'
  []
  [stress2]
    type = StokesLinearViscous
    mu = 2
    block = '2 3 5 8'
  []
[]

[BCs]
  [Periodic]
    [p]
      variable = p
      auto_direction = 'x y z'
    []
    [u]
      variable = u
      auto_direction = 'x y z'
    []
  []
[]

[Postprocessors]
  [avg_stress_xx]
     type = ElementAverageValue
     variable = stress_xx
  []
  [avg_stress_yy]
     type = ElementAverageValue
     variable = stress_yy
  []
  [avg_stress_zz]
     type = ElementAverageValue
     variable = stress_zz
  []
  [avg_stress_xy]
     type = ElementAverageValue
     variable = stress_xy
  []
  [avg_stress_yz]
     type = ElementAverageValue
     variable = stress_yz
  []
  [avg_stress_xz]
     type = ElementAverageValue
     variable = stress_xz
  []
[]

[Executioner]
  type = Steady

  solve_type = 'newton'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_mat_solver_package'
  petsc_options_value = 'lu NONZERO superlu_dist'
  line_search = none

  l_max_its = 10
  l_tol = 1e-8
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
[]
