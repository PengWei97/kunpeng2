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
#include "GrainDataTracker.h"

// Forward Declarations
class EulerAngleProvider;

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class PWComputePolycrystalElasticityTensor : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();
  // 共有数据接口

  PWComputePolycrystalElasticityTensor(const InputParameters & parameters);
  // 构造函数

protected:
  virtual void computeQpElasticityTensor();
  // 定义一个函数

  std::string _h_name;

  std::string _D_h_name;

  std::string _x_name;

  Real _length_scale;
  Real _pressure_scale;

  /// Grain tracker object
  const GrainDataTracker<RankFourTensor> & _grain_tracker;
  // 获取grain_tracker

  /// Number of order parameters
  const unsigned int _op_num;
  // 序参数数目

  /// Order parameters
  const std::vector<const VariableValue *> _vals;
  // 序参数的变量

  /// vector of elasticity tensor material properties
  std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;
  // 存放RankFourTensor类型的指针

  std::vector<MaterialProperty<Real> *> _h;

  std::vector<MaterialProperty<Real> *> _D_h;

  std::vector<MaterialProperty<Real> *> _x;

  /// Conversion factor from J to eV
  const Real _JtoeV;
};
