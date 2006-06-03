/* $Id: BCJKineticEqn.h,v 1.5 2004-07-15 08:29:14 paklein Exp $ */
#ifndef _BCJ_KINETIC_EQN_H_
#define _BCJ_KINETIC_EQN_H_

#include "KineticEqnBase.h"
#include "dArrayT.h"

namespace Tahoe {

class EVPFDBaseT;

/** kinetic equations for the BCJ model.
 * Implementation of the kinetic equations for the BCJ model.
 * Kinetic Equation:
   \f[
		\dot{\epsilon}_{eqp} = f(\sigma_{eff}, \kappa) = f \sinh \left[ \frac{\sigma_{eff} - \kappa - Y}{V} \right] 
   \f]
 * Dynamic Yield Condition:
   \f[
		\sigma_{eff} = h(\dot{\epsilon}_{eqp}, \kappa)
   \f]
 * Material parameters:
   \f[
		V = C_1 \exp \left( -\frac{C_2}{\Theta} \right)
   \f]
   \f[
		Y = \frac{C_3}{2 \left[ C_{21} + \exp \left(-\frac{C_4}{\Theta} \right) \right]} \left[1 + \tanh \left( C_{19}(C_{20} - \Theta) \right) \right]   
   \f]
   \f[
		f = C_5 \exp \left(-\frac{C_6}{\Theta} \right)
   \f]
 * with material constants: 
   \f[
		C_1 , C_2, C_3, C_4, C_5, C_6, C_{19}, C_{20}, C_{21}
   \f]
 * and temperature \f$ \Theta \f$.
 */
class BCJKineticEqn : public KineticEqnBase
{
 public:

  /** constructor */
  BCJKineticEqn(EVPFDBaseT& model);

  /** destructor */
  ~BCJKineticEqn();

  /** static yield stress */
  virtual double g(double eqp);

  /** \name kinetic equation functions */
  /*@{*/
  virtual double f        (double sigma, double kappa);
  virtual double DfDsigma (double sigma, double kappa);
  virtual double DfDs     (double sigma, double kappa);
  /*@}*/

  /** \name dynamic yield condition functions */
  /*@{*/
  virtual double h         (double eqpdot, double kappa);
  virtual double DhDeqpdot (double eqpdot, double kappa);
  virtual double DhDs      (double eqpdot, double kappa);
  /*@}*/

 private:

  /** material properties for Kinetic Eqn */
  void ComputeMaterialProperties(const double theta);

 private:

  /** temperature */
  double fTheta;

  /** \name material constants */
  /*@{*/
  double fC1, fC2, fC3, fC4, fC5, fC6;
  double fC19, fC20, fC21;
  /*@}*/
};

} // namespace Tahoe 
#endif  /* _BCJ_KINETIC_EQN_H_ */
