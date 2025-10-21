"""multi_padlock_design public API.

Exposes legacy `lib` modules at top level for backward compatibility while
preserving the ability to `import lib` directly.
"""

from importlib import import_module as _import_module

# Re-export everything from lib for existing notebooks that did
# `from multi_padlock_design import *` (since old root __init__ did from lib import *)
try:
    _lib = _import_module("lib")
    for _name in getattr(_lib, "__all__", []):
        globals()[_name] = getattr(_lib, _name)
    # If __all__ not defined, optionally expose common modules
    if not getattr(_lib, "__all__", None):
        for _name in dir(_lib):
            if not _name.startswith("_"):
                globals()[_name] = getattr(_lib, _name)
except ModuleNotFoundError:
    pass

__all__ = [n for n in globals() if not n.startswith("_")]
