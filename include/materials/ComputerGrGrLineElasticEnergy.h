#pragma once

#include "DerivativeFunctionMaterialBase.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"
#include "GrainDataTracker.h"

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class ComputerGrGrLineElasticEnergy : public DerivativeFunctionMaterialBase
{
public:
  static InputParameters validParams();

  ComputerGrGrLineElasticEnergy(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  // virtual void computerQpGrGrElasticityEnergy();
  virtual Real computeF() override;
  // virtual Real computeDF(unsigned int op_num) override;

  const std::string _base_name;

  std::string _elasticity_energy_name;

  std::string _elastic_energy_name;

  std::string _h_name;

  std::string _D_h_name;

  const MaterialProperty<RankTwoTensor> & _elastic_strain;
  const MaterialProperty<RankTwoTensor> & _stress;
  // const MaterialProperty<RankTwoTensor> & _lag_e_grgr;

  Real _length_scale;
  Real _pressure_scale;

  const GrainDataTracker<RankFourTensor> & _grain_tracker;

  /// Number of order parameters
  const unsigned int _op_num;

  /// Order parameters
  const std::vector<const VariableValue *> _vals;

  MaterialProperty<Real> & _elasticity_energy;

  MaterialProperty<Real> & _elastic_energy;

  std::vector<MaterialProperty<Real> *> _D_elastic_energy;

  std::vector<MaterialProperty<Real> *> _h;

  std::vector<MaterialProperty<Real> *> _D_h;

  /// Conversion factor from J to eV
  const Real _JtoeV;
};

