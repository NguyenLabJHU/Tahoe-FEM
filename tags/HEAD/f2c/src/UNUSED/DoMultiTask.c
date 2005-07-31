
/****************************************************************************
*																			*
*	This code implements a simple but flexible function	that can be			*
*	called to add multitasking to Fortran code ported to the Macintosh		*
*	using Mac F2C and the MetroWerks, THINK, or	Symantec C/C++ compilers.	*
*																			*
*	This code is based extensively on code provided by:						*
*																			*
*		Mel	Martinez														*
*		The	Johns Hopkins University										*
*		Dept. of Physics													*
*		mem@jhu.edu															*
*																			*
*	However, any and all bugs present are my responsiblity.					* 
*		IMT	16Sep95															*
*																			*
****************************************************************************/






/***************************************************************************

	Metrowerks CodeWarrior version
	******************************

	The approach is simple.  At program startup, F2Cmain() calls InitMultiTask()
	which set up a counter that keeps track of the ticks elapsed.   It also 
	initializes the spinning cursors (if the required resource is present; 
	otherwise it fails without impacting anything). 

	Then whenever you call the domultitask_ (or DoMultiTask) function, it checks
	if enough ticks have elapsed.  If so it calls WNE and spins the cursor.  If
	WNE returns an event, it is passed to the SIOUX console to handle.

	At program completion, F2Cmain() calls EndMultiTask() which cleans up.

***************************************************************************/
	
#ifdef CW_F2C_MAC

#include <Events.h>
#include <OSUtils.h>
#include <Quickdraw.h>
#include <ToolUtils.h>
#include <Resources.h>

#include <SIOUX.h>   


/* Declare f2c library function used to recover from exception exits & aborts */
void sig_die( char*, int );

/* Cursor related functions */
Boolean InitAnimatedCursors( short acurID );
void ReleaseAnimatedCursors( void );
void SpinMyCursor( void );

/* Private globals */
static long			gNextCheck;			/* Used to hold tick counts */
static long			gTickSlice = 2;		/* 1/30th second in ticks (each tick = 1/60th sec) */


/**********************************************************************
*																	  *
*	The	slice parameter above is the control that sets exactly 	      *
*	how often the program will actually	call WaitNextEvent.			  *
*																	  *
*		- The above	value is the recommended rate for assuring		  *
*		  that user	interaction	is not affected	(even QuickTime		  *
*		  movies should	still run fine in the foreground).			  *
*																	  *
*		- If you want to be	less friendly and use more CPU			  *
*		  time,	increase this number.								  *
*																	  *
**********************************************************************/




/*
	Use this in F2Cmain.c to initialize the multi-tasking code
*/

void InitMultiTask( long sliceInTicks ) 
{
	if ( sliceInTicks > 0 )
		gTickSlice = sliceInTicks;	

	InitAnimatedCursors( 128 );
	
	/* Start the tick count */
	gNextCheck = TickCount() + gTickSlice;
}





/*
	Use this in F2Cmain.c to close-out the multi-tasking code
*/

void EndMultiTask( void )
{
	ReleaseAnimatedCursors();
}





/*
	Get to this call in your FORTRAN code by inserting CALL DOMULTITASK( sleepTime )
*/

int domultitask_( long *sleepTime )
{
	extern Boolean SIOUXQuitting;
	EventRecord myEvent;

 	if( TickCount() > gNextCheck )   					/* Time to check for events again? */
    {
     	if( WaitNextEvent( everyEvent, &myEvent, *sleepTime, NULL ) )
      	{
      		/* Restore arrow while we're handling a real event */
      		SetCursor( &(qd.arrow) );
      		
      		/* Need to do something with the event if we got one */
			/* Add additional event handling code here if you need it */
			SIOUXHandleOneEvent( &myEvent );      
			if( SIOUXQuitting )
    			sig_die("User interrupt; execution stopped", 1);		/* Graceful quit */
 		}
    	SpinMyCursor(); 								/* Spin the cursor */
		gNextCheck = TickCount() + gTickSlice;			/* Reset the tick count */
    }
    
    return 1;
}



/* 
	Add this function so users can add DoMultiTask(n) calls to the
	C output from Mac F2C if they prefer, instead of inserting
	CALL DOMULTITASK(n) in the FORTRAN code 
*/

void DoMultiTask( long sleepTime )
{
	domultitask_( &sleepTime );
}






/********************

	Following are functions related to making the cursors spin.

	Wrote my own instead of using the InitCursors()/SpinCursor()/RotateCursor()
	package because it is an incomplete port from the MPW environment found in
	the MPW ToolLibs libraries and linking against them in the incorrect order
	can give rise to subtle inconsistencies.  A likely source of problems for users.

********************/



typedef struct		/* My version of the structure of an 'acur' resource */
{
	short numberOfFrames;			/* Number of cursors to animate */
	short whichFrame;				/* Current frame number	*/
	CursHandle frame[];				/* Pointer to the first cursor */
} acur, *acurPtr, **acurHandle;

static acurHandle gCursorList;		/* The cursor list */



/* 
	Try to get the acur record and the cursor list for acurID, 
	returning 1 (= TRUE) if everything goes as planned. 
*/

Boolean InitAnimatedCursors( short acurID )
{
	register short i = 0;
	register short cursID;
	Boolean noErrFlag = 0;

	if ( gCursorList = (acurHandle) GetResource( 'acur',acurID ) ) 
	{
		/* Got it! */ 
		noErrFlag = 1;
		while ( (i < (*gCursorList)->numberOfFrames) && noErrFlag ) 
		{
			/* The id of the cursor is stored in the high word of the frame handle */ 
			cursID = (short) HiWord( (long) (*gCursorList)->frame[i] );
			(*gCursorList)->frame[i] = GetCursor( cursID );
			if ( (*gCursorList)->frame[i] )
				i++;								/* Get the next one */
			else
				noErrFlag = 0;						/* Couldn't find the cursor */
		}
	}
	
	if ( noErrFlag ) 
	{
		(*gCursorList)->whichFrame = 0;
	}
	else
	{
		/* Free up memory we won't use */
		if ( gCursorList )
		{
			(*gCursorList)->numberOfFrames = i;		/* These are all we managed to get */
			ReleaseAnimatedCursors();
		}
		gCursorList = NULL;
	}
	
	return noErrFlag;
}




/* 
	Free up the storage used by the current animated cursor 
	and all of its frames 
*/

void ReleaseAnimatedCursors( void )
{
	short i;
	if ( gCursorList )
	{
		for ( i = 0; i< (*gCursorList)->numberOfFrames; i++ )
			ReleaseResource( (Handle) (*gCursorList)->frame[i] );
		ReleaseResource( (Handle) gCursorList );
	}
}



/* 	
	Display the next frame in the sequence. 
*/

void SpinMyCursor( void )
{
	if ( gCursorList )
	{
		/* Grab the frame, increment (and reset, if necessary) the count, */
		/* and display the new cursor */
		SetCursor( *((*gCursorList)->frame[(*gCursorList)->whichFrame++]) );
		if( (*gCursorList)->whichFrame == (*gCursorList)->numberOfFrames )
			(*gCursorList)->whichFrame = 0;
	}
}


#endif	/* CW_F2C_MAC */






/***************************************************************************

	THINK (68K) and Symantec (PPC) version
	**************************************

	The THINK console does not have to pass an event back to it (i.e., an 
	equivalent of CW's SIOUXHandleOneEvent function).  So we can't call
	WNE directly because if we get an event there' no way to hand it back 
	to the console (sorry, setting the event mask in WNE to 0x0000 or 0x0001
	to get only null events doesn't work--nice idea though).  

	So instead have to do the following:
	(a) Set stdin to raw, unbuffered mode.
	(b) Call getchar() which forces the console to call WNE for us.  
        Because we are in raw mode, getchar() returns immediately even if
        the user doesn't type anything.
	(c) Set stdin back to buffered mode.

	Alas, changing the mode of stdin is very slow.  So instead of implementing
	the above directly, we make raw/unbuffer the 'default' mode on stdin so 
	we don't have to change modes just to call WNE (via getchar).  And we
	compensate by changing back to buffered mode prior to doing any input on 
	stdio (and then restoring raw/unbuffered mode).  (InitMultiTask() and
	EndMultiTask() are used to adjust the default mode on stdio).

	Luckily for this plan, the f2c library bottlenecks all input operations
	through fread(), getc(), and ungetc().  Our replacements for these are
	defined below.  They simply pass the call through to the 'real' versions 
	of the above functions, changing modes appropriately if the file used 
	for I/O is stdin.

	My versions of these functions are patched into the f2c library code
	via macros in TPM_I.h and TPM_F.h.

	Many thanks to Phil Shapiro of Symantec for suggesting this approach.

	Also, there is no cursor spinning in this TPM/SPM version because the 
	console updates the cursor on every opportunity, so all we do is end 
	up fighting the console for control of the cursor's appearance.

***************************************************************************/
	

#if defined(TPM_F2C) || defined(SPM_F2C)

/* Undefine some #defines that come from TPM_I.h and SPM_I.h that we don't want here */
#undef fread
#undef __getc
#undef ungetc
#include <stdio.h>			/* We want the original unsubstituted versions */

#include <Events.h>
#include <OSEvents.h>
#include <OSUtils.h>

#include <console.h>

/* Private globals */
static long			gNextCheck = 0;		/* Used to hold tick counts */
static long			gTickSlice = 2;		/* 1/30th second in ticks (each tick = 1/60th sec) */

/**********************************************************************
*																	  *
*	The	slice parameter above is the control that sets exactly 	      *
*	how often the program will actually	call WaitNextEvent.			  *
*																	  *
*		- The above	value is the recommended rate for assuring		  *
*		  that user	interaction	is not affected	(even QuickTime		  *
*		  movies should	still run fine in the foreground).			  *
*																	  *
*		- If you want to be	less friendly and use more CPU			  *
*		  time,	increase this number.								  *
*																	  *
**********************************************************************/




/*
	This is used in F2Cmain.c to initialize the multi-tasking code
*/

void InitMultiTask( long sliceInTicks ) 
{
	if ( sliceInTicks > 0 )
		gTickSlice = sliceInTicks;	
	
	/* Set console so we can interogate for input (to trigger WNE) without waiting for input */
	csetmode( C_RAW, stdin );

	/* Start the tick count */
	gNextCheck = TickCount() + gTickSlice;
}





/*
	This is used in F2Cmain.c to close-out the multi-tasking code
*/

void EndMultiTask( void )
{
    csetmode(C_ECHO, stdin); 	/* Restore standard, buffered mode (not really required :) */
}





/*
	Get to this call in your FORTRAN code by inserting CALL DOMULTITASK( sleepTime )
*/

int domultitask_( long *sleepTime )
{
 	if( TickCount() > gNextCheck )   				/* Time to check for events again? */
    {
    	getchar();   								/* poll for char to have console call WNE */
		gNextCheck = TickCount() + gTickSlice;		/* Reset the tick count */
    }
    return 1;
}



/* 
	Add this function so users can add DoMultiTask(n) calls to the
	C output from Mac F2C if they prefer, instead of inserting
	CALL DOMULTITASK(n) in the FORTRAN code 
*/

void DoMultiTask( long sleepTime )
{
	domultitask_( &sleepTime );
}



/*
	Because now stdin is (by default) in raw mode, must re-set it before reading
    from it.  These replacement functions cover all the input functions used 
	in the f2c libs.  #defines in TPM_F2C and SPM_F2C ensure that these are called
    instead of their originals
*/


size_t F2C_fread( void *ptr, size_t size, size_t n, FILE *f )
{
	size_t  nr;
	if ( f == stdin )
	{
		csetmode(C_ECHO, stdin); 
		nr = fread( ptr, size, n, f );
		csetmode(C_RAW, stdin);
	}
	else
		nr = fread( ptr, size, n, f );

	return  nr;
}


int F2C_getc( FILE * f)
{
	int i;
	if ( f == stdin )
	{
		csetmode(C_ECHO, stdin); 
		i = __getc( f );
		csetmode(C_RAW, stdin);
	}
	else
		i = __getc( f );

	return  i;
}


int F2C_ungetc( int c, FILE *f )
{
	int i;
	if ( f == stdin )
	{
		csetmode(C_ECHO, stdin); 
		i = ungetc( c, f );
		csetmode(C_RAW, stdin);
	}
	else
		i = ungetc( c, f );

	return  i;
}
	

#endif	/* TPM_F2C or SPM_F2C */



