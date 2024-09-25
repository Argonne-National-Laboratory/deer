# Simple 3D test

[GlobalParams]
  displacements = "disp_x disp_y disp_z"
[]

[Mesh]
  type = FileMesh
  file = "small.exo"
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]
[UserObjects]
  [euler_angle_file]
    type = ElementPropertyReadFile
    nprop = 3
    prop_file_name = 2gr_Euler.tex
    read_type = block
    nblock = 2
    blocks_zero_numbered = false
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
    boundary = back
    variable = disp_y
    value = 0.0
  []
  [leftz]
    type = DirichletBC
    preset = true
    boundary = bottom
    variable = disp_z
    value = 0.0
  []
  [pull_z]
    type = FunctionDirichletBC
    boundary = top
    variable = disp_z
    function = pfn
    preset = true
  []
[]

[Functions]
  [pfn]
    type = PiecewiseLinear
    x = '0    10'
    y = '0.00 0.1'
  []
[]

[AuxVariables]
  [orientation_q1]
    order = CONSTANT
    family = MONOMIAL
  []
  [orientation_q2]
    order = CONSTANT
    family = MONOMIAL
  []
  [orientation_q3]
    order = CONSTANT
    family = MONOMIAL
  []
  [orientation_q4]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [q1]
    type = MaterialStdVectorAux
    property = orientation
    index = 0
    variable = orientation_q1
  []
  [q2]
    type = MaterialStdVectorAux
    property = orientation
    index = 1
    variable = orientation_q2
  []
  [q3]
    type = MaterialStdVectorAux
    property = orientation
    index = 2
    variable = orientation_q3
  []
  [q4]
    type = MaterialStdVectorAux
    property = orientation
    index = 3
    variable = orientation_q4
  []
[]

[Materials]
  [stress1]
    type = NEMLCrystalPlasticity
    database = "test.xml"
    model = "grain_1"
    large_kinematics = true
    euler_angle_reader = euler_angle_file
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        add_variables = true
        new_system = true
        formulation = TOTAL
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
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

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  start_time = 0
  dt = 1
  end_time = 10.0
[]

[Outputs]
  exodus = true
[]
