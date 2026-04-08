//------------------------------------------------------------------------------
// CHOLMOD/Utility/t_cholmod_error: CHOLMOD error handling
//------------------------------------------------------------------------------

// CHOLMOD/Utility Module. Copyright (C) 2023, Timothy A. Davis, All Rights
// Reserved.
// SPDX-License-Identifier: LGPL-2.1+

//------------------------------------------------------------------------------

#include "cholmod_internal.h"

//------------------------------------------------------------------------------
// cholmod_error
//------------------------------------------------------------------------------

int CHOLMOD(error)
(
    // input:
    int status,             // Common->status
    const char *file,       // source file where error occurred
    int line,               // line number where error occurred
    const char *message,    // error message to print
    cholmod_common *Common
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    RETURN_IF_NULL_COMMON (FALSE) ;

    //--------------------------------------------------------------------------
    // set the error status
    //--------------------------------------------------------------------------

    Common->status = status ;

    //--------------------------------------------------------------------------
    // handle the error, unless we're inside a CHOLMOD try/catch block
    //--------------------------------------------------------------------------

    if (!(Common->try_catch))
    {

        //----------------------------------------------------------------------
        // call the user error handler, if present
        //----------------------------------------------------------------------

        if (Common->error_handler != NULL)
        {
            Common->error_handler (status, file, line, message) ;
        }
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    return (TRUE) ;
}
