//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTrackerAdd.h"
#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"

class ComputerGrGrElasticEnergy : public DerivativeMaterialInterface<Material>
{
  public: 
    static InputParameters validParams();

    ComputerGrGrElasticEnergy(const InputParameters & parameters);

  protected:
    virtual void computerQpGrGrElasticityEnergy();

    const std::string _base_name;
    std::string _elasticity_energy_name;

    const MaterialProperty<RankTwoTensor> & _pk2_grgr;
    const MaterialProperty<RankTwoTensor> & _lag_e_grgr;
    MaterialProperty<RankTwoTensor> & _piaolak2;

    Real _length_scale;
    Real _pressure_scale;

    const GrainDataTrackerAdd<RankFourTensor,RealVectorValue> & _grain_tracker;

    /// Number of order parameters
    const unsigned int _op_num;

    /// Order parameters
    const std::vector<const VariableValue *> _vals;

    MaterialProperty<Real> & _elasticity_energy;

    // MaterialProperty<RankTwoTensor> & _piaola2;


    std::vector<MaterialProperty<Real> *> _D_elastic_energy;

    /// Conversion factor from J to eV
    const Real _JtoeV;
};