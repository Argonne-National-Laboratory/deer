# Test under stress cotrolled condition BC. Load reversal is very quick

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -1
  zmax = 1
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 0'
    top_right = '0.5 0.5 0.5'
  []
  [./scale]
  type = TransformGenerator
  input = new_block
  transform = SCALE
  vector_value ='0.06 0.06 0.06'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = scale
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./applied_disp_z]
    type = PiecewiseLinear
    x = '0.0 1000 1100 1500'
    y = '0 0.06 0 -0.03'
  [../]
[]

[BCs]
  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./z1]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function  = applied_disp_z
  [../]
[]

# Constraint System
[Constraints]
  [./x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'right'    # boundary
    penalty = 1e6
  [../]
  [./y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'top'    # boundary
    penalty = 1e6
  [../]
[]

[NEMLMechanics]
  displacements = 'disp_x disp_y disp_z'
  kinematics = large
  add_all_output = true
  add_displacements = true
  formulation = total
[]


[Modules/TensorMechanics/CohesiveZoneMaster]
  [czm]
    strain = FINITE
    boundary = 'interface'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
  []
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "mat.xml"
    model = "creep_and_hardening"
    large_kinematics = true
  [../]
  [./czm_mat]
    type = GBCavitation
    boundary = 'interface'
    max_time_cut = 10
    D_failure = 0.9
    max_nonlinear_iter = 20
    minimum_allowed_stiffness = 1
    minimum_allowed_residual_life = 10
    nucleation_on = true
    growth_on = true
    use_old_bulk_property = false
    nl_residual_abs_tol = 1e-12
    interface_thickness = 1e-4
    E_penalty_after_failure_minus_thickenss = 1e8
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  # Executioner
  type = Transient
  automatic_scaling = true
  compute_scaling_once = false
  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5
  l_max_its = 2
  nl_max_its = 10
  start_time = 0.0
  [./TimeStepper]
   type =IterationAdaptiveDT
    optimal_iterations = 10
    dt =  100
  []
  dtmin = 1e-10
  end_time = 1500
  dtmax = 50
[]

[Outputs]
  exodus = true
  print_nonlinear_converged_reason = false
  print_linear_converged_reason = false
[]
