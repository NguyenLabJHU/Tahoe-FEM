/* $Id: ArgSpecT.h,v 1.4 2002-07-05 22:26:24 paklein Exp $ */

#ifndef _ARG_SPEC_T_H_
#define _ARG_SPEC_T_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** definition of command arguments */
class ArgSpecT
{
  public:
  
 	/** type specifications */
	enum ArgTypeT {int_ = 0, 
	        double_ = 1, 
	        string_ = 2, 
	          bool_ = 3, 
             float_ = 4};

	/** constructor for unnamed argument */
	ArgSpecT(ArgTypeT);

	/** constructor for named argument */
	ArgSpecT(ArgTypeT type, const StringT& name);

	/** copy constructor */
	ArgSpecT(const ArgSpecT& arg);

	/** destructor */
	~ArgSpecT(void);

	/** return the argument name */
	const StringT& Name(void) const { return fName; };

	/** return the argument type */
	const ArgTypeT& Type(void) const { return fType; };
	
	/** return the type name */
	const char* TypeName(void) const;

	/** return the argument name */
	const StringT& Prompt(void) const { return fPrompt; };

	/** set argument prompt */
	void SetPrompt(const StringT& prompt) { fPrompt = prompt; };
	
	/** return true argument has default value */
	bool HasDefault(void) const { return fDefault != NULL; };
	
	/** clear the default value */
	void ClearDefault(void);
	
	/* set default value */
	void SetDefault(const int&);
	void SetDefault(const double&);
	void SetDefault(const StringT&);
	void SetDefault(const bool&);
	void SetDefault(const float&);
	void SetDefault(const char* s) { SetDefault(StringT(s)); };

	/* set default value */
	void GetDefault(int&) const;
	void GetDefault(double&) const;
	void GetDefault(StringT&) const;
	void GetDefault(bool&) const;
	void GetDefault(float&) const;

	/** return true argument has been assigned a value */
	bool HasValue(void) const { return fValue != NULL; };

	/** clear the default value */
	void ClearValue(void);
	
	/** read value from stream.
	 * \return true if successful, false if read fails */
	bool ReadValue(istream& in);

	/* set default value */
	void SetValue(const int&);
	void SetValue(const double&);
	void SetValue(const StringT&);
	void SetValue(const bool&);
	void SetValue(const float&);
	void SetValue(const char* s) { SetValue(StringT(s)); };

	/* set default value */
	void GetValue(int&) const;
	void GetValue(double&) const;
	void GetValue(StringT&) const;
	void GetValue(bool&) const;
	void GetValue(float&) const;
	
	/** set the value to its default. Has not effect if the value
	 * has not default. */
	void SetValueToDefault(void);

	/** assignment operator */
	ArgSpecT& operator=(const ArgSpecT& rhs);

	/** write spec to output */
	void Write(ostream& out) const;

	/** write value to the output stream */
	void WriteValue(ostream& out) const;

  private:
  
  	/** command name */
  	StringT fName;
  	
  	/** argument type */
  	ArgTypeT fType;  	

	/** prompt */
	StringT fPrompt;
	
	/** default value */
	void* fDefault;

	/** value */
	void* fValue;
};

} // namespace Tahoe 
#endif /* _ARG_SPEC_T_H_ */
