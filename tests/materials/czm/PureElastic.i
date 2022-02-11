[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  zmax = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 1'
    top_right = '1 1 2'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = false
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy '
                          'cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy '
                          'mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
  [../]
[]
[Materials]
  [./stress]
    type = CauchyStressFromNEML
    database = "../../neml_test_material.xml"
    model = "elastic_model"
  [../]
  [./czm]
    type = PureElasticCZM
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    E = 1e2
    G = 1e2
    interface_thickness = 1
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
  nl_abs_tol = 1e-6
  dtmin = 1
  dtmax = 1
  end_time = 6
[]

[BCs]
  [./x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
  [./x_top]
    type = DirichletBC
    boundary = front
    variable = disp_x
    value = 0.0
  [../]
  [./y_top]
    type = DirichletBC
    boundary = front
    variable = disp_y
    value = 0.0
  [../]
  [./z_top]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = disp_fun
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 1  3  4 6'
    y = '0 1 -1 -2 0'
  [../]
[]



[Outputs]
  exodus = true
  sync_times = '1 2 4 6'
[]
