// $Id: NOX_Tahoe_Vector.h,v 1.5 2002-07-05 22:28:43 paklein Exp $
#ifndef NOX_TAHOE_VECTOR_H
#define NOX_TAHOE_VECTOR_H

/* optional */
#ifdef __NOX__

// base class
#include "NOX_Abstract_Vector.H"

namespace Tahoe {

// forward declarations
class dArrayT;

//! %NOX %Tahoe support.

//! Implementation of NOX::Abstract::Vector for %Tahoe vectors.
class Vector : public NOX::Abstract::Vector {

  public:			

	//! Default constructor. Makes empty vector
	Vector(void);

	//! Construct by copying map and/or elements of an dArrayT.
	Vector(const dArrayT& source, NOX::CopyType type = NOX::DeepCopy);

	//! Destruct Vector.
	~Vector();

	//@{ \name Access to underlying toolbox vector.

	//! type conversion to underlying toolbox vector.
	operator dArrayT&() { return *fArray; };
	
	//! type conversion to underlying toolbox vector.
	operator const dArrayT&() const { return *fArray; };

	//! Get reference to underlying toolbox vector.
	dArrayT& get_dArrayT(void) { return *fArray; };
	
	//! Get const reference to underlying toolbox vector.
	const dArrayT& get_dArrayT() const { return *fArray; };

  //@}

  //@{ \name Initialization methods.

  virtual NOX::Abstract::Vector& init(double value);

  //! assignment operator
  virtual NOX::Abstract::Vector& operator=(double value) { return init(value); } ;

  //! Copies source vector into "this".
  virtual NOX::Abstract::Vector& operator=(const dArrayT& source);

  virtual NOX::Abstract::Vector& operator=(const Vector& source);
  //! See above.
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& source);
  
  virtual NOX::Abstract::Vector& abs(const Vector& source);
  //! See above.
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& source);

  virtual NOX::Abstract::Vector& reciprocal(const Vector& source);
  //! See above.
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& source);

  //@}

  //@{ \name Update methods.

  virtual NOX::Abstract::Vector& scale(double gamma);

  virtual NOX::Abstract::Vector& update(double alpha, const Vector& a, 
			     double gamma = 0.0);
  //! See above.
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double gamma = 0.0);

  virtual NOX::Abstract::Vector& update(double alpha, const Vector& a, 
			     double beta, const Vector& b,
			     double gamma = 0.0);
  //! See above.
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double beta, const NOX::Abstract::Vector& b,
			     double gamma = 0.0);

  //@}

  //@{ \name Creating new Vectors. 

  virtual NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;

  //@}

  //@{ \name Norms.

  virtual double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TWO) const;

  virtual double norm(const Vector& weights, NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TWO) const;
  //! See above.
  virtual double norm(const NOX::Abstract::Vector& weights, 
		      NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TWO) const;

  //@}

  //@{ \name Dot products

  virtual double dot(const Vector& y) const;
  //! See above.
  virtual double dot(const NOX::Abstract::Vector& y) const;

  //@}

  virtual int length() const;

	/** dimension the vector based on the source. Copies the shape of the source. */
	void DimensionTo(const dArrayT& source);

 private:
  
  //! Pointer to dArrayT owned by this object
  dArrayT* fArray;

};
} // namespace Tahoe 
#endif /* __NOX__ */
#endif /* NOX_TAHOE_VECTOR_H */
