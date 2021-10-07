# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 00:07:52 2021

@author: veenstra
"""

import os
import pandas as pd
from inspect import getmembers, isfunction
import importlib
import hatyan

print('## generating copy contents of README.md to docstring of hatyan/__init__')
dir_scripts = os.path.dirname(os.path.realpath(__file__))
file_init = os.path.join(dir_scripts,'..','hatyan/__init__.py')
file_readme = os.path.join(dir_scripts,'..','README.md')

with open(file_readme) as f:
    data_readme = pd.Series(f.readlines())
init_pre = pd.Series(['# -*- coding: utf-8 -*-\n',
                      '# this __init__.py file is automatically generated with scripts/generate_documentation.py\n',
                      '"""\n'])
init_post = pd.Series(['\n',
                       '"""\n',
                       '\n',
                       '__author__ = """Jelmer Veenstra"""\n',
                       "__email__ = 'jelmer.veenstra@deltares.nl'\n","__version__ = '%s'\n"%(hatyan.__version__),
                       '\n',
                       'from hatyan.analysis_prediction import *\n',
                       'from hatyan.astrog import *\n',
                       'from hatyan.components import *\n',
                       'from hatyan.foreman_core import *\n',
                       'from hatyan.hatyan_core import *\n',
                       'from hatyan.timeseries import *\n',
                       'from hatyan.wrapper_RWS import *\n'])

data_init_out = pd.concat([init_pre,data_readme,init_post])
data_init_out_str = data_init_out.str.cat()

with open(file_init,'w') as f:
    f.write(data_init_out_str)

print('## generating HTML documentation')
err_html = os.system('pdoc --html hatyan -o doc --force --config sort_identifiers=False')

#print('## generating pdf documentation')
#err_pdf_pdoc = os.system('pdoc --pdf hatyan --config sort_identifiers=False > doc/hatyan_functions.md')
#err_pdf_pandoc = os.system('pandoc --metadata=title:"hatyan-%s documentation, automatically generated from script comments and function docstrings" --toc --toc-depth=4 --variable fontsize=12pt --from=markdown+abbreviations --pdf-engine=xelatex --top-level-division=chapter --output=doc/hatyan_functions.pdf doc/hatyan_functions.md'%(hatyan.__version__))

if err_html: # or err_pdf_pdoc or err_pdf_pandoc:
    raise Exception('ERROR: in documentation generation, check feedback above')
else:
    print('## succesfully finished')
