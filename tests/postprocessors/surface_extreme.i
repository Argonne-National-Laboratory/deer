[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        new_system = true
        formulation = UPDATED
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[Materials]
  [stress]
    type = CauchyStressFromNEML
    database = "../test_materials.xml"
    model = "elastic_model"
    large_kinematics = false
  []
  [tensor_rate]
    type = TensorRateMaterial
    rank_two_tensor = cauchy_stress
  []
[]

[BCs]
  [x_0]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [y_0]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [z_0]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  [z_1]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
    value = 0.1
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [max_pull]
    type = SideExtremePostprocessor
    variable = disp_z
    boundary = 'front'
    value_type = max
  []

  [min_pull]
    type = SideExtremePostprocessor
    variable = disp_z
    boundary = 'front'
    value_type = min
  []

  [max_right]
    type = SideExtremePostprocessor
    variable = disp_x
    boundary = 'right'
    value_type = max
  []

  [min_right]
    type = SideExtremePostprocessor
    variable = disp_x
    boundary = 'right'
    value_type = min
  []
[]

[Executioner]
  type = Transient

  end_time = 1.0
  dt = 1.0

  solve_type = NEWTON
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 4
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]
[Outputs]
  csv = true
[]
