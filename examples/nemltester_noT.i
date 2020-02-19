[Mesh]
      file = 'one.exo'
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
      active = 'TensorMechanics'
      [./TensorMechanics]
            displacements = 'disp_x disp_y disp_z'
      [../]
[]

[Materials]
      [./strain]
            type = ComputeSmallStrain
            displacements = 'disp_x disp_y disp_z'
      [../]
      [./stress]
            type = ComputeNEMLStress
            database = 'materials.xml'
            model = 'chaboche_600'
      [../]
[]

[Functions]
      [./straining]
            type = PiecewiseLinear
            x = '0.0 1.0'
            y = '0.0 0.1'
      [../]
[]

[BCs]
      [./left]
            type = DirichletBC
            preset = true
            variable = disp_x
            value = 0.0
            boundary = 'left'
      [../]
      [./back]
            type = DirichletBC
            preset = true
            variable = disp_y
            value = 0.0
            boundary = 'back'
      [../]
      [./bottom]
            type = DirichletBC
            preset = true
            variable = disp_z
            value = 0.0
            boundary = 'bottom'
      [../]
      [./top]
            type = FunctionDirichletBC
            variable = disp_z
            boundary = 'top'
            function = straining
            preset = true
      [../]
[]

[Executioner]
      type = Transient
      solve_type = PJFNK

      petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
      petsc_options_value = 'hypre boomeramg 101'

      nl_max_its = 20
      nl_rel_tol = 1.0e-6
      nl_abs_tol = 1.0e-10
      l_max_its = 200
      l_tol = 1e-10

      start_time = 0.0
      end_time = 1.0

      dt = 1.0e-2
[]

[AuxVariables]
      [./s11]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./s22]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./s33]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./s12]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./s13]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./s23]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e11]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e22]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e33]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e12]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e13]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./e23]
            order = CONSTANT
            family = MONOMIAL
      [../]

      [./von_mises]
            order = CONSTANT
            family = MONOMIAL
      [../]
[../]

[AuxKernels]
      [./report_s11]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 0
            index_j = 0
            variable = s11
      [../]

      [./report_s22]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 1
            index_j = 1
            variable = s22
      [../]

      [./report_s33]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 2
            index_j = 2
            variable = s33
      [../]

      [./report_s12]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 0
            index_j = 1
            variable = s12
      [../]

      [./report_s13]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 0
            index_j = 2
            variable = s13
      [../]

      [./report_s23]
            type = RankTwoAux
            rank_two_tensor = stress
            index_i = 1
            index_j = 2
            variable = s23
      [../]

      [./report_e11]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 0
            index_j = 0
            variable = e11
      [../]

      [./report_e22]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 1
            index_j = 1
            variable = e22
      [../]

      [./report_e33]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 2
            index_j = 2
            variable = e33
      [../]

      [./report_e12]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 0
            index_j = 1
            variable = e12
      [../]

      [./report_e13]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 0
            index_j = 2
            variable = e13
      [../]

      [./report_e23]
            type = RankTwoAux
            rank_two_tensor = total_strain
            index_i = 1
            index_j = 2
            variable = e23
      [../]

      [./von_mises_kernel]
            type = RankTwoScalarAux
            variable = von_mises
            rank_two_tensor = stress
            scalar_type = VonMisesStress
      [../]
[../]

[Outputs]
      exodus = true
      print_perf_log = true
[]
