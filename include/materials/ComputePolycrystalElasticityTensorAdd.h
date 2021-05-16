//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeElasticityTensorBase.h"
#include "GrainDataTrackerAdd.h"

// Forward Declarations
class EulerAngleProvider;

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class ComputePolycrystalElasticityTensorAdd : public ComputeElasticityTensorBase 
{
public:
  static InputParameters validParams();

  ComputePolycrystalElasticityTensorAdd(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  Real _length_scale;
  Real _pressure_scale;

  /// Grain tracker object
  const GrainDataTrackerAdd<RankFourTensor,RealVectorValue> & _grain_tracker;

  /// Number of order parameters
  const unsigned int _op_num;

  /// Order parameters
  const std::vector<const VariableValue *> _vals;

  /// vector of elasticity tensor material properties
  std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;

  /// Crystal Rotation Matrix
  MaterialProperty<RankTwoTensor> & _crysrot;

  /// Conversion factor from J to eV
  const Real _JtoeV;
};
