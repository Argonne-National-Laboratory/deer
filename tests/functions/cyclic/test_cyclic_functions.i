[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    nx = 1
    ny = 1
    nz = 1
    dim = 3
  [../]
  [./node0_1]
    input = msh
    type = ExtraNodesetGenerator
    new_boundary = 'node01'
    nodes = '0 1'
  [../]
  [./node0]
    input = node0_1
    type = ExtraNodesetGenerator
    new_boundary = 'node0'
    nodes = '0'
  [../]
  [./node67]
    input = node0
    type = ExtraNodesetGenerator
    new_boundary = 'node67'
    nodes = '6 7'
  [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    family = LAGRANGE
    order = FIRST
  []
  [./disp_y]
    family = LAGRANGE
    order = FIRST
  []
  [./disp_z]
    family = LAGRANGE
    order = FIRST
  []
[]

[Kernels]
  [./sdx]
      type = StressDivergenceNEML
      variable = disp_x
      component = 0
      block = 0
  [../]
  [./sdy]
      type = StressDivergenceNEML
      variable = disp_y
      component = 1
      block = 0
  [../]
  [./sdz]
      type = StressDivergenceNEML
      variable = disp_z
      component = 2
      block = 0
  [../]
[]

[Functions]
  [./cycle_time_fun]
    type = CycleTime
    cycle_period = 3
  [../]
  [./cyclic_displacment_fun]
    type = PiecewiseLinearCycle
    x = '0 1 3'
    y = '0 1.5e-3 0'
    cycle_time_func = cycle_time_fun
  [../]
  [./cycle_fraction_fun]
    type = CycleFraction
    cycle_period = 3
  [../]
  [./cycle_numeber_fun]
    type = CycleNumber
    cycle_period = 3
  [../]
[]

[AuxVariables]
  [./cycle_time]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cycle_fraction]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cycle_number]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cyclic_displacment]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./cycle_time_aux]
    type = FunctionAux
    function = cycle_time_fun
    variable = cycle_time
  [../]
  [./cyclic_displacment_aux]
    type = FunctionAux
    function = cyclic_displacment_fun
    variable = cyclic_displacment
  [../]
  [./cycle_fraction_aux]
    type = FunctionAux
    function = cycle_fraction_fun
    variable = cycle_fraction
  [../]
  [./cycle_numeber_aux]
    type = FunctionAux
    function = cycle_numeber_fun
    variable = cycle_number
  [../]
[]

[BCs]
  [./x0_x]
    type = DirichletBC
    variable = disp_x
    boundary = node0
    value = 0.0
  [../]
  [./x0_y]
    type = DirichletBC
    variable = disp_y
    boundary = node01
    value = 0.0
  [../]
  [./x0_z]
    type = DirichletBC
    variable = disp_z
    boundary = node01
    value = 0.0
  [../]


  [./z1_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = node67
    function = cyclic_displacment_fun
  [../]
  [./z1_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = node67
    function = cyclic_displacment_fun
  [../]
[]





[Materials]
  [./strain]
    type = ComputeNEMLStrain
  [../]
  [./stress1]
    type = ComputeNEMLStressUpdate
    database = "simple_elastic.xml"
    model = "simple_elastic"
    block = 0
  [../]
[]


 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
   [../]
 []


[Executioner]

    # Preconditisoned JFNK (default)
    type = Transient
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'lu superlu_dist'

    solve_type = NEWTON
    nl_abs_tol = 1e-8
    nl_rel_tol = 1e-8
    nl_max_its = 7
    l_tol = 1e-15
    # l_max_its = 5
    start_time = 0.0
    dtmin = 1e-50
    dtmax = 1
    end_time = 9

  line_search = none
[]
[Outputs]
  sync_times = '0.5 1'
  exodus =true
[]
