#ifndef _BAGON_ASSERT_H_
#define _BAGON_ASSERT_H_

#ifdef MEX_COMPILE

#include "mex.h"
#define assert(x) if (!(x)) { mexErrMsgTxt(#x);}

#else

#define assert(x) if (!(x)) { printf(#x); exit(1);}

#endif



#endif // _BAGON_ASSERT_H_
