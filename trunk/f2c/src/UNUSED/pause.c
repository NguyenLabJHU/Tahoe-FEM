/*
	Create a MacOS pause function.
	
	IMT  14Aug97
*/


#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC)

#include <Events.h>

int pause( void )
{
	while ( !Button() )
	{
		SystemTask();
	}
	
	return 1;
}

#endif /* MacOS compilers */
