[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  zmax = 2
  elem_type = HEX20
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 0 1'
    top_right = '1 1 2'
  []
  [./boundary]
    type = SideSetsBetweenSubdomainsGenerator
    input = new_block
    master_block = '0 1'
    paired_block = 1
    new_boundary = 'interface'
  []
  [node]
    type = ExtraNodesetGenerator
    coord = '0.0 0.0 0.0'
    new_boundary = 6
    input = boundary
  []
  [snode]
    type = ExtraNodesetGenerator
    coord = '1.0 0.0 0.0'
    new_boundary = 7
    input = node
  []
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
 [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]



[Materials]
  [./fplastic]
    type = FiniteStrainPlasticMaterial
    block = '0 1'
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '0 1'
    C_ijkl = '2.827e5 1.21e5 1.21e5 2.827e5 1.21e5 2.827e5 0.808e5 0.808e5 0.808e5'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = '0 1'
    displacements = 'disp_x disp_y disp_z'
  [../]

  [./vonmises_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = vonmises_interface
    effective_stress_type = VonMises
    boundary = 'interface'
  [../]
  [./s1_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = s1_interface
    effective_stress_type = MaxPrincipal
    boundary = 'interface'
  [../]
  [./tresca_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = tresca_interface
    effective_stress_type = Tresca
    boundary = 'interface'
  [../]
  [./hydrostatic_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = hydrostatic_interface
    effective_stress_type = Hydrostatic
    boundary = 'interface'
  [../]
  [./huddleston_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = huddleston_interface
    params_vector = '0.03'
    effective_stress_type = Huddleston
    boundary = 'interface'
  [../]
  [./hayhurst_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = hayhurst_interface
    params_vector = '0.2 0.3'
    effective_stress_type = Hayhurst
    boundary = 'interface'
  [../]
  [./rccmrx_mises_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = rccmrx_mises_interface
    params_vector = '0.2'
    effective_stress_type = RCCMRXMises
    boundary = 'interface'
  [../]
  [./rccmrx_tresca_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = rccmrx_tresca_interface
    params_vector = '0.2'
    effective_stress_type = RCCMRXTresca
    boundary = 'interface'
  [../]
  [./max_s1_mises_interface]
    type = EffectiveStressInterfaceMaterial
    effective_stress_mp_name = max_s1_mises_interface
    effective_stress_type = maxS1AndMises
    boundary = 'interface'
  [../]
[]

[Functions]
 [./topfunc]
   type = ParsedFunction
   value = '0.01*t'
 [../]
[]

[BCs]
  [./bottom3]
    type = DirichletBC
    variable = disp_z
    boundary = 0
    value = 0.0
    preset = true
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 5
    function = topfunc
    preset = true
  [../]
  [./corner1]
    type = DirichletBC
    variable = disp_x
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./corner2]
    type = DirichletBC
    variable = disp_y
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./corner3]
    type = DirichletBC
    variable = disp_z
    boundary = 6
    value = 0.0
    preset = true
  [../]
  [./side1]
    type = DirichletBC
    variable = disp_y
    boundary = 7
    value = 0.0
    preset = true
  [../]
  [./side2]
    type = DirichletBC
    variable = disp_z
    boundary = 7
    value = 0.0
    preset = true
  [../]
[]

[AuxVariables]
  [./vonmises_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s1_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tresca_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hydrostatic_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./huddleston_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hayhurst_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rccmrx_mises_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rccmrx_tresca_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./max_s1_mises_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]

  [./vonmises_interface]
    type = MaterialRealAux
    variable = vonmises_interface
    property = vonmises_interface
    boundary = 'interface'
  [../]
  [./s1_interface]
    type = MaterialRealAux
    variable = s1_interface
    property = s1_interface
    boundary = 'interface'
  [../]
  [./tresca_interface]
    type = MaterialRealAux
    variable = tresca_interface
    property = tresca_interface
    boundary = 'interface'
  [../]
  [./hydrostatic_interface]
    type = MaterialRealAux
    variable = hydrostatic_interface
    property = hydrostatic_interface
    boundary = 'interface'
  [../]
  [./huddleston_interface]
    type = MaterialRealAux
    variable = huddleston_interface
    property = huddleston_interface
    boundary = 'interface'
  [../]
  [./hayhurst_interface]
    type = MaterialRealAux
    variable = hayhurst_interface
    property = hayhurst_interface
    boundary = 'interface'
  [../]
  [./rccmrx_mises_interface_interface]
    type = MaterialRealAux
    variable = rccmrx_mises_interface
    property = rccmrx_mises_interface
    boundary = 'interface'
  [../]
  [./rccmrx_tresca_interface_interface]
    type = MaterialRealAux
    variable = rccmrx_tresca_interface
    property = rccmrx_tresca_interface
    boundary = 'interface'
  [../]
  [./max_s1_mises_interface_interface]
    type = MaterialRealAux
    variable = max_s1_mises_interface
    property = max_s1_mises_interface
    boundary = 'interface'
  [../]
[]

[Postprocessors]

  [./vonmises_interface]
    type = SideAverageValue
    variable = vonmises_interface
    boundary = 'interface'
  [../]
  [./s1_interface]
    type = SideAverageValue
    variable = s1_interface
    boundary = 'interface'
  [../]
  [./tresca_interface]
    type = SideAverageValue
    variable = tresca_interface
    boundary = 'interface'
  [../]
  [./hydrostatic_interface]
    type = SideAverageValue
    variable = hydrostatic_interface
    boundary = 'interface'
  [../]
  [./huddleston_interface]
    type = SideAverageValue
    variable = huddleston_interface
    boundary = 'interface'
  [../]
  [./hayhurst_interface]
    type = SideAverageValue
    variable = hayhurst_interface
    boundary = 'interface'
  [../]
  [./rccmrx_mises_interface]
    type = SideAverageValue
    variable = rccmrx_mises_interface
    boundary = 'interface'
  [../]
  [./rccmrx_tresca_interface]
    type = SideAverageValue
    variable = rccmrx_tresca_interface
    boundary = 'interface'
  [../]
  [./max_s1_mises_interface]
    type = SideAverageValue
    variable = max_s1_mises_interface
    boundary = 'interface'
  [../]
[]

[Executioner]

  type = Transient

  dt=0.1
  dtmin=0.1
  dtmax=1
  end_time=1.0

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
