## @package pysub
# Python Submission scripts
#
# In principle running the EUTelescope processors to get your data
# analyzed is not so difficult, but one has to remember which
# processors are best suited for a particular procedure, the correct
# order and reasonable parameter default values.
#
# Moreover if you think that Marlin and EUTelescope can also be
# executed on the GRID this is introducing an extra degree of
# complexity!
#
# To simplify all this, we are now providing together with EUTelescope
# itselft and the eutelHistogram package, also a set of very powerful
# submission script, written in python with the main goal to allow any
# user or even a sofware robot to submit his/her own job either
# locally or on the GRID without having to modify a single line of
# steering file!
#
# Of course this come at a price! Marlin and EUTelescope were born to
# be as flexibile as possible but all this flexibility is making an
# automatic job submission nearly impossible.
#
# pysub is a sort of compromise between automatization and
# flexibility. The idea is that all the main parameters of pysub
# itself, like the GRID configuration or the verbosity level can be
# adjusted changing their value in a simple ASCII configuration file
# (config/config.cfg). The execution parameters instead are given in
# some steering file used as templates by the script.
#
# Our advice is the following:
#  @li Everytime you start a new analysis
#  copy the config/config.cfg into a new file, say
#  config/MyConfig.cfg.
#  @li In the new configuration properly set all the GRID paramters in
#  the GRID section.
#  @li Again in the configuration file modify the section concerning the
#  steering template file name in order to point to your own modified
#  version, if this is required.
#  @li Tell pysub to use your own config file and not the default one
#  either specifing the configuration file as a command line option
#  (--config-file = CONFIGFILE ) or by setting the enviromental
#  variable SUBMIT_CONFIG.
#
# From a programmer point of view, the pysub package is designed like
# this:
# <h4> The base class </h4>
#  The base submitter class (SubmitBase) is a very simple class
# contaning just the overall scheme of all the operation that a
# submitter has to have. Those are <b>configure</b>, <b>execute</b>
# and <b>end</b>.
#
# <h4> The submitter class </h4>
#  This is a class inheriting from SubmitBase, like SubmitConverter or
# SubmitTest implementing the call back the base class.
#
# <h4> The submitter script </h4>
#  By it self the submitter class is not running, but just defining a
# functionality that has to be executed. This is the main job of the
# submitter script like submit-converter.py that is actually
# responsible to create a new instance of the submitter class and make
# it working.
#
# @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
# @version $Id: __init__.py,v 1.13 2009-06-04 17:16:50 bulgheroni Exp $
#
from submitbase              import SubmitBase
from submitconverter         import SubmitConverter
from submittest              import SubmitTest
from submitnativecopy        import SubmitNativeCopy
from submitpedestal          import SubmitPedestal
from submitclusearch_nopede  import SubmitCluSearch
from submitfilter_nopede     import SubmitFilter
from submiteta_nopede        import SubmitEta
from submithitmaker_nopede   import SubmitHitMaker
from submitalign_nopede      import SubmitAlign
from submitfitter_nopede     import SubmitFitter
from submitanadut            import SubmitAnaDUT
from submitapplyAlignment    import SubmitApplyAlignment
from error                   import *

