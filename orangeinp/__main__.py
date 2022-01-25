# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
"""
Generate an org.xml from an input org.omn file.
"""
from .run import run

get_parser = run.make_argparser
main = run.main

if __name__ == '__main__':
    main()
