#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
# --------------------------
# custom package
# --------------------------

from PATTYpipe.Utility      import (sp,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   createDIR)

# --------------------------
# main 
# --------------------------
def stepFinal_summary(conf_dict,logfile):
    wlog('Collect results',logfile)
    tmpresult = 'tmpResults/'
    createDIR(tmpresult)
    sp("mv tmp_*.* %s"%(tmpresult))

    if conf_dict['General']['mode'] == "bulk":    
        sp("mv %s_correctSig.bdg %s"%(conf_dict['General']['outname'],tmpresult))
    else:
        sp("mv tmp_scProcess/ %s"%(tmpresult))
    if conf_dict['options']['keeptmp']:
        wlog('--keeptmp was set, keep intermediate results',logfile)
        pass
    else:
        wlog('--keeptmp was not set, remove intermediate results',logfile)
        sp("rm -r tmpResults/")

    return conf_dict

















