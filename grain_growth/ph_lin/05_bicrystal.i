# 03: 弹性能打包传输到kernel


[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 3
  xmax = 20
  ymax = 20
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = 10
      y2 = 20
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  # [./elastic_strain11]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./elastic_strain22]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./elastic_strain12]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./C1111]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./active_bounds_elemental]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./euler_angle]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  # [./elastic_strain11]
  #   type = RankTwoAux
  #   variable = elastic_strain11
  #   rank_two_tensor = elastic_strain
  #   index_i = 0
  #   index_j = 0
  #   execute_on = timestep_end
  # [../]
  # [./elastic_strain22]
  #   type = RankTwoAux
  #   variable = elastic_strain22
  #   rank_two_tensor = elastic_strain
  #   index_i = 1
  #   index_j = 1
  #   execute_on = timestep_end
  # [../]
  # [./elastic_strain12]
  #   type = RankTwoAux
  #   variable = elastic_strain12
  #   rank_two_tensor = elastic_strain
  #   index_i = 0
  #   index_j = 1
  #   execute_on = timestep_end
  # [../]
  # [./unique_grains]
  #   type = FeatureFloodCountAux
  #   variable = unique_grains
  #   flood_counter = grain_tracker
  #   execute_on = 'initial timestep_begin'
  #   field_display = UNIQUE_REGION
  # [../]
  # [./var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   flood_counter = grain_tracker
  #   execute_on = 'initial timestep_begin'
  #   field_display = VARIABLE_COLORING
  # [../]
  # [./C1111]
  #   type = RankFourAux
  #   variable = C1111
  #   rank_four_tensor = elasticity_tensor
  #   index_l = 0
  #   index_j = 0
  #   index_k = 0
  #   index_i = 0
  #   execute_on = timestep_end
  # [../]
  # [./active_bounds_elemental]
  #   type = FeatureFloodCountAux
  #   variable = active_bounds_elemental
  #   field_display = ACTIVE_BOUNDS
  #   execute_on = 'initial timestep_begin'
  #   flood_counter = grain_tracker
  # [../]
  # [./euler_angle]
  #   type = OutputEulerAngles
  #   variable = euler_angle
  #   euler_angle_provider = euler_angle_file
  #   grain_tracker = grain_tracker
  #   output_euler_angle = 'phi1'
  # [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    # value = if(t<=19.0,-0.01*t,0.19)
    value = 0.1
  [../]
[]

[BCs]
  # [./top_displacement]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = top
  #   value = -10.0
  # [../]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = tdisp
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 4 # nm
#     GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    GBMobility = 1.0e-12
    time_scale = 1.0e-9
    length_scale = 1.0e-9
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
  [./elasticenergy] 
    type = GetMaterialParamsLine
    # args = 'gr0 gr1' 
    grain_tracker = grain_tracker
    outputs = exodus 
  [../]

[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = test.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    connecting_threshold = 0.05
    compute_var_to_feature_map = true
    flood_entity_type = elemental
    execute_on = 'initial timestep_begin'

    euler_angle_provider = euler_angle_file
    fill_method = symmetric9
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'

    outputs = none
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0,gr1 disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

  start_time = 0.0
  end_time = 20
  # num_steps = 300
  dt = 0.2

  [./Adaptivity]
   initial_adaptivity = 2
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 2
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
