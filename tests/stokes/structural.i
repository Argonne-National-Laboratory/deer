[Mesh]
  second_order = true
  [mesh]
    type = CartesianMeshGenerator
    dim = 3
    dx = '1 1'
    dy = '1 1'
    dz = '2'
    ix = '1 1'
    iy = '1 1'
    iz = '2'
    subdomain_id = '1 2 3 4'
  []
[]

[Variables]
  [p]
    order = FIRST
    family = MONOMIAL
  []
  [u]
    order = SECOND
    family = LAGRANGE_VEC
  []
[]

[AuxVariables]
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
[]

[AuxKernels]
    [stress_xx]
        type = ADMaterialRankTwoTensorAux
        variable = stress_xx
        property = 'stress'
        i = 0
        j = 0
    []
    [stress_yy]
        type = ADMaterialRankTwoTensorAux
        variable = stress_yy
        property = 'stress'
        i = 1
        j = 1
    []
    [stress_zz]
        type = ADMaterialRankTwoTensorAux
        variable = stress_zz
        property = 'stress'
        i = 2
        j = 2
    []
[]

[ICs]
  [u]
    type = VectorConstantIC
    variable = u
    x_value = 1e-15
    y_value = 1e-15
    z_value = 1e-15
  []
[]

[Kernels]
  [equil]
    type = ADStokesStressDivergence
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
    mu = 1.0
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

  [test]
    type = ADVectorFunctionNeumannBC
    boundary = 'front'
    variable = u
    function_x = '0'
    function_y = '0'
    function_z = '2.0/3.0'
  []
[]

[Executioner]
  type = Steady

  solve_type = 'newton'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
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
