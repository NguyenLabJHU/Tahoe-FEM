// $Id: NOX_Tahoe_Vector.h,v 1.2 2002-04-02 23:30:55 paklein Exp $
#ifndef NOX_TAHOE_VECTOR_H
#define NOX_TAHOE_VECTOR_H

/* optional */
#ifdef __NOX__

// base class
#include "NOX_Abstract_Vector.H"

// forward declarations
class dArrayT;

namespace NOX {

//! %NOX %Tahoe support.
namespace Tahoe {

//! Implementation of NOX::Abstract::Vector for %Tahoe vectors.
class Vector : public Abstract::Vector {

  public:			

	//! Default constructor. Makes empty vector
	Vector(void);

	//! Construct by copying map and/or elements of an dArrayT.
	Vector(const dArrayT& source, CopyType type = DeepCopy);

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

  virtual Abstract::Vector& init(double value);

  //! assignment operator
  virtual Abstract::Vector& operator=(double value) { return init(value); } ;

  //! Copies source vector into "this".
  virtual Abstract::Vector& operator=(const dArrayT& source);

  virtual Abstract::Vector& operator=(const Vector& source);
  //! See above.
  virtual Abstract::Vector& operator=(const Abstract::Vector& source);
  
  virtual Abstract::Vector& abs(const Vector& source);
  //! See above.
  virtual Abstract::Vector& abs(const Abstract::Vector& source);

  virtual Abstract::Vector& reciprocal(const Vector& source);
  //! See above.
  virtual Abstract::Vector& reciprocal(const Abstract::Vector& source);

  //@}

  //@{ \name Update methods.

  virtual Abstract::Vector& scale(double gamma);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double gamma = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double gamma = 0.0);

  virtual Abstract::Vector& update(double alpha, const Vector& a, 
			     double beta, const Vector& b,
			     double gamma = 0.0);
  //! See above.
  virtual Abstract::Vector& update(double alpha, const Abstract::Vector& a, 
			     double beta, const Abstract::Vector& b,
			     double gamma = 0.0);

  //@}

  //@{ \name Creating new Vectors. 

  virtual Abstract::Vector* clone(CopyType type = DeepCopy) const;

  //@}

  //@{ \name Norms.

  virtual double norm(NormType type = TWO) const;

  virtual double norm(const Vector& weights, NormType type = TWO) const;
  //! See above.
  virtual double norm(const Abstract::Vector& weights, NormType type = TWO) const;

  //@}

  //@{ \name Dot products

  virtual double dot(const Vector& y) const;
  //! See above.
  virtual double dot(const Abstract::Vector& y) const;

  //@}

  virtual int length() const;

	/** dimension the vector based on the source. Copies the shape of the source. */
	void DimensionTo(const dArrayT& source);

 private:
  
  //! Pointer to dArrayT owned by this object
  dArrayT* fArray;

};
} // namespace Tahoe
} // namespace NOX

#endif /* __NOX__ */
#endif /* NOX_TAHOE_VECTOR_H */
