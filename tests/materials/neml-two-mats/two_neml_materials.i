[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  []
  [./subdomain]
    input = msh
    type = ElementSubdomainIDGenerator
    subdomain_ids = '0 0 0 0
                     1 1 1 1'
  [../]
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
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy '
                          'cauchy_stress_xz cauchy_stress_yz strain_xx strain_yy strain_zz strain_xy '
                          'strain_xz strain_yz'
      []
    []
  []
[]

[BCs]
  [./back_x]
    type = DirichletBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
  [./back_y]
    type = DirichletBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]
  [./back_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./front_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = 1e-3*t
  [../]
[]


[Materials]
  [./stress_0]
    type = CauchyStressFromNEML
    database = neml_model.xml
    model = model_1
    block = 0
  [../]
  [./stress_1]
    type = CauchyStressFromNEML
    database = neml_model.xml
    model = model_2
    block = 1
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  solve_type = NEWTON

  dt = 0.5
  end_time = 1
  line_search = none
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
