[Mesh]
  second_order = true
  [mesh]
    type = CartesianMeshGenerator
    dim = 3
    dx = '1 1'
    dy = '1 1'
    dz = '2'
    ix = '2 2'
    iy = '2 2'
    iz = '4'
    subdomain_id = '1 2 3 4'
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
  [stress]
    type = StokesLinearViscous
    mu = 1
  []
[]

[BCs]
  [fix_x]
    type = ADVectorFunctionDirichletBC
    variable = u
    boundary = 'left'
    set_y_comp = False
    set_z_comp = False
  []

  [fix_y]
    type = ADVectorFunctionDirichletBC
    variable = u
    boundary = 'bottom'
    set_x_comp = False
    set_z_comp = False
  []

  [fix_z]
    type = ADVectorFunctionDirichletBC
    variable = u
    boundary = 'back'
    set_x_comp = False
    set_y_comp = False
  []

  [fix_p1]
    type = ADDirichletBC
    variable = p
    boundary = 'front'
    value = -1.0
  []

  [fix_zero]
    type = ADDirichletBC
    variable = p
    boundary = 'right top'
    value = 0.0
  []
[]

[Executioner]
  type = Steady

  solve_type = 'newton'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
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
