/*
	This function returns 0 if file called fileName exists,
	1 otherwise.
	
	IMT 28 Nov 94
*/


#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC)

#include <Files.h>

static StringPtr CopyCtoPstr( StringPtr pStr, char *cStr );

int f2c_access( char *fileName, int notUsed )
{	
	OSErr 	err;
	short 	vRefNum;
	FInfo 	info;
	Str255	pFileName;
	
	CopyCtoPstr( pFileName, fileName );			/* Make name a Pascal string */
		
	err = GetFInfo( pFileName, 0, &info );		/* Check on the default volume */
	
	if ( err )
		return  1;
	else
		return  0;
}


static StringPtr CopyCtoPstr( StringPtr pStr, char *cStr )
{
	short		i;
	char		*p;
	
	i = 0;
	p = ((char *) pStr) + 1;
	while ( *cStr && i < 255 )
	{
		*p++ = *cStr++;
		i++;
	}
	*pStr = i;
	
	return  (StringPtr) pStr;
}


#endif	/* MacOS compilers */

