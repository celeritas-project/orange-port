# -*- coding: utf-8 -*-
# Copyright 2021 UT-Battelle, LLC and SCALE Developers.
# See the top-level COPYRIGHT file for details.
from omnutils.db.script import Script, Regenerator

run = Script('orange', "ORANGE model", 'orangeinp', suffix='.org')
generate_xml = Regenerator(run)
