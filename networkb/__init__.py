import networkb.classes
from networkb.classes import *

import networkb.algorithms
from networkb.algorithms import *

import networkb.ploters
from networkb.ploters import *

import networkb.IO
from networkb.IO import *

import logging
logger = logging.getLogger('networkb')
#if len(logger.handlers) == 0:	# To ensure reload() doesn't add another one
#    logger.addHandler(NullHandler())

try:
  import pkg_resources
  __version__ = pkg_resources.get_distribution("networkb").version
except:
    __version__ = '?'
