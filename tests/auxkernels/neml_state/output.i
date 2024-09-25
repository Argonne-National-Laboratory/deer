# Simple 3D test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = TOTAL
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[AuxVariables]
  [alpha]
    order = CONSTANT
    family = MONOMIAL
  []
  [X1_0]
    order = CONSTANT
    family = MONOMIAL
  []
  [X2_0]
    order = CONSTANT
    family = MONOMIAL
  []
  [X3_0]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [alpha]
    type = NEMLStateAux
    database = "../../test_materials.xml"
    model = "chaboche_model"
    variable = alpha
    state_variable = "alpha"
  []
  [X1_0]
    type = NEMLStateAux
    database = "../../test_materials.xml"
    model = "chaboche_model"
    variable = X1_0
    state_variable = "backstress_0_0"
  []
  [X2_0]
    type = NEMLStateAux
    database = "../../test_materials.xml"
    model = "chaboche_model"
    variable = X2_0
    state_variable = "backstress_1_0"
  []
  [X3_0]
    type = NEMLStateAux
    database = "../../test_materials.xml"
    model = "chaboche_model"
    variable = X3_0
    state_variable = "backstress_2_0"
  []
[]

[Functions]
  [pullz]
    type = ParsedFunction
    expression = '0.005 * t'
  []
[]

[BCs]
  [leftx]
    type = DirichletBC
    preset = true
    boundary = left
    variable = disp_x
    value = 0.0
  []
  [lefty]
    type = DirichletBC
    preset = true
    boundary = bottom
    variable = disp_y
    value = 0.0
  []
  [leftz]
    type = DirichletBC
    preset = true
    boundary = back
    variable = disp_z
    value = 0.0
  []
  [pull_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = pullz
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "../../test_materials.xml"
    model = "chaboche_model"
    large_kinematics = false
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 0.1
  dtmin = 0.1
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
