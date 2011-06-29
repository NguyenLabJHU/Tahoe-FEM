/* $Id: NOX_Tahoe_Group.h,v 1.6 2002-07-06 02:13:23 paklein Exp $ */
#ifndef NOX_TAHOE_GROUP_H
#define NOX_TAHOE_GROUP_H

/* optional */
#ifdef __NOX__

/* base class */
#include "NOX_Abstract_Group.H"

/* direct members */
#include "NOX_Tahoe_Vector.h" // for NOX::CopyType
#include "NOX_Common.H" // class data element (string)
#include "dArrayT.h"

/* forward declarations */
namespace NOX {
	namespace Parameter {
		class List;
	}
}
namespace Tahoe {
  class SolverInterfaceT;
  class GlobalMatrixT;
}

namespace Tahoe {

//! %Tahoe support. 

/** implementation of NOX::Abstract::Group for %Tahoe groups */
class Group: public NOX::Abstract::Group {

public:
  
	/** constructor. The resulting object uses, but does not take ownership
	 * of the matrix passed in. */
	Group(SolverInterfaceT& interface, dArrayT& X, GlobalMatrixT& J);

	/** copy constructor. The resulting object owns the Jacobian matrix used 
	 * by the class */
	Group(const Group& source, NOX::CopyType type);

	/** destructor */
	virtual ~Group(void);
	
	/** access to the solution vector 
  /*! 
    \brief Copies the values of all vectors and any other data in
    source group to this group.  (May invalidate shared data for
    source group.)
  */
  virtual NOX::Abstract::Group& operator=(const Group& source);
      
  //! See above.
  virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      
   /** \name "Compute" functions */
  //@{ 

  //! Compute and return solution vector, x, where this.x = grp.x() + step * d.
  virtual bool computeX(const Group& grp, const Vector& d, double step);

  //! See above.
  virtual bool computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);

  //! Compute and return RHS. Also computes and stores norm of RHS.
  virtual bool computeRHS(void);

  //! Compute Jacobian.
  virtual bool computeJacobian(void);

  //! Compute and return gradient.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual bool computeGrad(void);

  //! Compute and return Newton direction, using parameters for nonlinear solve.
  //! Throws an error if RHS and Jacobian have not been computed.
  virtual bool computeNewton(NOX::Parameter::List& params);

  //@}

  /** @name Jacobian operations.
   
    Operations using the Jacobian matrix. These may not be defined in
    matrix-free scenarios.
  */

  //@{
  
  /*! 
    \brief If supported, returns true and calculates 
    result = Jacobian * input.  Otherwise, returns false.  
    Returns false is any errors occur, such as the Jacobian not being
    computed.
  */
  virtual bool applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

  /*! 
    \brief If supported, returns true and calculates result =
    Jacobian^T * input.  Otherwise, returns false.  Throws an error if
    the Jacobian has not been computed. 
  */
  virtual bool applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  
  /*!
    \brief Applies the Jacobian Diagonal to the given input vector.
  */
  virtual bool applyJacobianDiagonalInverse(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
    

  //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns true
    if the corresponding "compute" function has been called since the
    last update to the solution vector (via instantiation or
    computeX).
  */

  //@{
  
  //! Return true if the RHS is valid.
  virtual bool isRHS() const { return fIsRHS; };
  //! Return true if the Jacobian is valid.
  virtual bool isJacobian() const { return fIsJacobian; };
  //! Return true if the gradient is valid.
  virtual bool isGrad() const { return fIsGrad; };
  //! Return true if the Newton direction is valid.
  virtual bool isNewton() const { return fIsNewton; };

  //@}

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const NOX::Abstract::Vector& getX() const { return fSolution; };

  //! Return right-hand-side (RHS). 
  virtual const NOX::Abstract::Vector& getRHS() const { return fRHS; };

  //! Return 2-norm of RHS.
  virtual double getNormRHS() const { return fRHSNorm; };

  //! Return gradient.
  virtual const NOX::Abstract::Vector& getGrad() const { return fGradient; };

  //! Return Newton direction.
  virtual const NOX::Abstract::Vector& getNewton() const { return fNewton; };

  //@}


  /** \name Creating new Groups. */
  //@{
  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  If type
    is "DeepCopy", then we need to create an exact replica of
    "this". Otherwise, if type is "CopyShape", we need only replicate
    the shape of "this". Returns NULL if clone is not supported.
  */
  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //@}

  private:
  
  	/** interface to Tahoe */
  	SolverInterfaceT& fTahoeInterface;
  
  	/** true if object owns its associated jacobian matrix */
  	bool fOwnJ;

	/** Jacobian matrix */
  	GlobalMatrixT* fJ;

	/** solution vector */
	Vector fSolution;

  	/** solution vector */
  	Vector fRHS;
	
	/** gradient of objective function */
	Vector fGradient;

	/** gradient of objective function */
	Vector fNewton;
	
	/* flags about what is current */
	bool fIsRHS;
	bool fIsJacobian;
	bool fIsGrad;
	bool fIsNewton;
	
	/** 2-norm of RHS. Computed during Group::computeRHS */
	double fRHSNorm;	
};

} // namespace Tahoe 
#endif /* __NOX__ */
#endif /* NOX_TAHOE_GROUP_H */
