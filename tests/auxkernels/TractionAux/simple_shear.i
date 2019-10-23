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
  [./loadFunction]
    type = ParsedFunction
    value = 1e-2*t
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
    function = loadFunction
  [../]
  [./z1_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = node67
    function = loadFunction
  [../]
[]

[AuxVariables]
  [./stress_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_zz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_yz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_xz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_xy]
    family = MONOMIAL
    order = CONSTANT
  [../]


  [./strain_xx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_zz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_yz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_xz]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./strain_xy]
    family = MONOMIAL
    order = CONSTANT
  [../]





  [./traction_N]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_S1]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./traction_S2]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./shear_traction_norm]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./traction_N]
    type = TractionAux
    scalar_type = 'normal'
    variable = traction_N
    property = stress
    boundary = 'front'
  [../]
  [./traction_S1]
    type = TractionAux
    scalar_type = 'shear1'
    variable = traction_S1
    property = stress
    boundary = 'front'
  [../]
  [./traction_S2]
    type = TractionAux
    scalar_type = 'shear2'
    variable = traction_S2
    property = stress
    boundary = 'front'
  [../]
  [./shear_traction_norm]
    type = TractionAux
    scalar_type = 'shear_norm'
    variable = shear_traction_norm
    property = stress
    boundary = 'front'
  [../]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 2
    index_j = 1
  []
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 2
    index_j = 0
  []
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  []


  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  []
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  []
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  []
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_yz
    index_i = 2
    index_j = 1
  []
  [./strain_xz]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xz
    index_i = 2
    index_j = 0
  []
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = mechanical_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  []
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
    end_time = 1

  line_search = none
[]
[Outputs]
  sync_times = '0.5 1'
  exodus =true
[]
