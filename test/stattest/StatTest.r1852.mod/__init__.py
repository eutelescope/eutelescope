__all__ = [ 'Interface','Tests','IO','ROOTIO','Utils' ]
__majorversion__ = 0
__minorversion__ = 1
__patchlevel__ = 0
__version__ = ".".join(str(e) for e in [__majorversion__,__minorversion__,__patchlevel__])

VERSION = __version__
def getMajorVersion(): return __majorversion__
def getMinorVersion(): return __minorversion__
def getPatchLevelVersion(): return __patchlevel__
def getVersion(): return __version__

