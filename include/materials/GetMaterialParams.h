#pragma once

#include "DerivativeFunctionMaterialBase.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"
#include "GrainDataTrackerAdd.h"

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class GetMaterialParams : public DerivativeFunctionMaterialBase
{
public:
  static InputParameters validParams();

  GetMaterialParams(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computerQpGrGrElasticityEnergy();
  virtual Real computeF() override;
  // virtual Real computeDF(unsigned int i_var) override;

  const std::string _base_name;

  std::string _elasticity_energy_name;

  const MaterialProperty<RankTwoTensor> & _pk2_grgr;
  const MaterialProperty<RankTwoTensor> & _lag_e_grgr;

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

  /// Stress tensor
  const MaterialProperty<RankTwoTensor> & _stress;

  ///@{ Strain and derivatives
  const MaterialProperty<RankTwoTensor> & _strain;
  ///@}
};