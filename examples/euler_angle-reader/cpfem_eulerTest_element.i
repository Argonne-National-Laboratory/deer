[Mesh]
  type = FileMesh
  file = 'n10.exo'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    order = second
  []
  [disp_y]
    order = second
  []
  [disp_z]
    order = second
  []
[]

[Functions]
  [pfn]
    type = PiecewiseLinear
    x = '0    100'
    y = '0.00 0.01'
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

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = FINITE
        new_system = true
        formulation = TOTAL
        volumetric_locking_correction = false
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_xy cauchy_stress_xz cauchy_stress_yz mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz'
      []
    []
  []
[]

[UserObjects]
  [euler_angle_file]
    type = ElementPropertyReadFile
    nprop = 3
    prop_file_name = grn_element_rand.tex
    read_type = element
  []
[]
[BCs]
  [left]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = left
    value = 0.0
  []

  [bottom]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bottom
    value = 0.0
  []

  [back]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = back
    value = 0.0
  []

  [front]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = pfn
    preset = true
  []
[]

[Materials]
  [stress]
    type = NEMLCrystalPlasticity
    database = "test.xml"
    model = "grain_1"
    large_kinematics = true
    euler_angle_reader = euler_angle_file
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

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_max_its = 20

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8

  end_time = 10.0
  dtmin = 0.1
  dt = 1.0
[]

[Outputs]
  exodus = true
[]
