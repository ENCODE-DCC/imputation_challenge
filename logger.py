#!/usr/bin/env python3
"""Imputation challenge logger

Author:
    Jin Lee (leepc12@gmail.com)
"""

import sys
import logging

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
