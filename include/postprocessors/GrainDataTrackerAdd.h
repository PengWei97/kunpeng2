//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainTracker.h"

/**
 * GrainTracker derived class template to base objects on which maintain physical
 * parameters for individual grains.
 */
template <typename T1,typename T2>
class GrainDataTrackerAdd : public GrainTracker
{
public:
  GrainDataTrackerAdd(const InputParameters & parameters);

  /// return data for selected grain
  const T1 & getDataElasticity(unsigned int grain_id) const;
  // 被computerElasticitytensor所调用，获取赋予给每个晶粒的弹性模量

  const T2 & getDataRotation(unsigned int grain_id) const;
  // 被computerElasticitytensor所调用，获取赋予给每个晶粒的旋转矩阵

protected:
  /// implement this method to initialize the data for the new grain
  virtual T1 newGrainElasticity(unsigned int new_grain_id) = 0;

  virtual T2 newGrainRotation(unsigned int new_grain_id) = 0;

  virtual void newGrainCreated(unsigned int new_grain_id);
  // 新晶粒的创建

  /// per grain data
  std::vector<T1> _grain_data_elasticity;

  std::vector<T2> _grain_data_rotation;
};

template <typename T1,typename T2>
GrainDataTrackerAdd<T1,T2>::GrainDataTrackerAdd(const InputParameters & parameters) : GrainTracker(parameters)
{
}

template <typename T1,typename T2>
const T1 &
GrainDataTrackerAdd<T1,T2>::getDataElasticity(unsigned int grain_id) const
{
  mooseAssert(grain_id < _grain_data_elasticity.size(), "Requested data for invalid grain index.");
  return _grain_data_elasticity[grain_id];
}

template <typename T1,typename T2>
const T2 &
GrainDataTrackerAdd<T1,T2>::getDataRotation(unsigned int grain_id) const
{
  mooseAssert(grain_id < _grain_data_elasticity.size(), "Requested data for invalid grain index.");
  return _grain_data_rotation[grain_id];
}

template <typename T1,typename T2>
void
GrainDataTrackerAdd<T1,T2>::newGrainCreated(unsigned int new_grain_id)
{
  if (_grain_data_elasticity.size() <= new_grain_id)
  {
    _grain_data_elasticity.resize(new_grain_id + 1);
    _grain_data_rotation.resize(new_grain_id + 1);

  }

  _grain_data_elasticity[new_grain_id] = newGrainElasticity(new_grain_id);

  _grain_data_rotation[new_grain_id] = newGrainRotation(new_grain_id);
}


