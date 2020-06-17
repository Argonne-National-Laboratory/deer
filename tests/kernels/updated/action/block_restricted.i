# Simple 3D test

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  []
  [./new_block]
     type = SubdomainBoundingBoxGenerator
     input = msh
     block_id = 1
     bottom_left = '0 0 0.5'
     top_right = '1 1 1'
  []
  [./add_boundary]
    type = BreakMeshByBlockGenerator
    input = new_block
  [../]
  [./lower_d_block]
    type = LowerDBlockFromSidesetGenerator
    input = add_boundary
    new_block_id = 10
    sidesets = 'interface'
  []
[]

[NEMLMechanics]
  displacements = "disp_x disp_y disp_z"
  kinematics = small
  add_all_output = true
  block = '0 1' #NEMLMechanics shoudl not be used on subdomain 10
[]

[Kernels] # required to add at least one kernel to block 10
  [null_x]
    type = NullKernel
    block = 10
    variable = 'disp_x'
  []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm1]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[BCs]
  [./leftx]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./lefty]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./leftz]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./pull_z]
    type = DirichletBC
    boundary = front
    variable = disp_z
    value = 0.1
  [../]
[]

[Materials]
  [./stress]
    type = ComputeNEMLStressUpdate
    database = "test.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]
  [./czm_3dc]
    type = SalehaniIrani3DCTraction
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 100
    maximum_shear_traction = 70
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
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

  dt = 1
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
