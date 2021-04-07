[Mesh]
  [msh]
    type = FileMeshGenerator
    # file = testfile.cpr
    has_fake_neighbors = true
    fake_neighbor_list_file_name = 'fake_neighbors_test_bmbb.csv'
  []
[]

[AuxVariables]
  [bmbb_element_id]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [set_material_id]
    type = ElemExtraIDAux
    variable = bmbb_element_id
    extra_id_name = bmbb_element_id
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]


[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]

[CohesiveZoneDeer]
   boundary = 'interface'
[]


[Functions]
  [applied_disp_z]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.1'
  []
[]


[BCs]
  [x0]
    type = DirichletBC
    variable = disp_x
    boundary = x0
    value = 0.0
  []
  [y0]
    type = DirichletBC
    variable = disp_y
    boundary = y0
    value = 0.0
  []
  [z0]
    type = DirichletBC
    variable = disp_z
    boundary = z0
    value = 0.0
  []
  [z1]
    type = FunctionDirichletBC
    boundary = z1
    function = applied_disp_z
    variable = disp_z
  []
[]

# Constraint System
[Constraints]
  [x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'x1'    # boundary
    penalty = 1e6
  []
  [y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'y1'    # boundary
    penalty = 1e6
  []
  [z1]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    secondary = 'z1'    # boundary
    penalty = 1e6
  []
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "../mat.xml"
    model = "elastic_model"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = PureElasticCZMIncremental
    E = 1e3
    G = 1e2
    boundary = 'interface'
  [../]
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -stol'
  petsc_options_value = 'lu superlu_dist 0'
  solve_type = NEWTON
  line_search = none
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 10
  l_tol = 1e-15
  l_max_its = 2
  dtmin = 0.1
  dtmax = 0.1
  end_time = 0.3
  dt = 0.1
[]


[Outputs]
  sync_times = '0 0.1 1 10 100 1000 10000 100000 1000000'
  csv = true
  #print_linear_converged_reason = false
  #print_nonlinear_converged_reason = false
[]
