[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 2
  ny = 4
  xmax = 3
  ymax = 2
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '1.5 0 0'
    top_right = '3 2 0'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
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
  [../]
[]

[Materials]
  [./stress]
    type = CauchyStressFromNEML
    database = "../../neml_test_material.xml"
    model = "elastic_model"
    large_kinematics = false
  [../]
  [./czm]
    type = PureElasticCZM
    boundary = 'interface'
    E = 1e0
    G = 1e0
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
  dt = 2
  end_time = 2
[]

[BCs]
  [./fix_x]
    type = DirichletBC
    boundary = left
    variable = disp_x
    value = 0.0
  [../]
  [./fix_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./move_x2]
    type = FunctionDirichletBC
    boundary = right
    variable = disp_x
    function = disp_fun
  [../]
[]

[Functions]
  [./disp_fun]
    type = PiecewiseLinear
    x = '0 2'
    y = '0 0.3'
  [../]
[]

[CZMStrain]
   boundary = interface
   strain = SMALL
   block = '0 1'
[]

[Outputs]
  csv = true
[]
