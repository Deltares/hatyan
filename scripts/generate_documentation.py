# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 00:07:52 2021

@author: veenstra
"""

import os

print('## generating HTML documentation')
err_html = os.system('pdoc --html hatyan -o doc --force --config sort_identifiers=False')

#print('## generating pdf documentation')
#err_pdf_pdoc = os.system('pdoc --pdf hatyan --config sort_identifiers=False > doc/hatyan_functions.md')
#err_pdf_pandoc = os.system('pandoc --metadata=title:"hatyan-%s documentation, automatically generated from script comments and function docstrings" --toc --toc-depth=4 --variable fontsize=12pt --from=markdown+abbreviations --pdf-engine=xelatex --top-level-division=chapter --output=doc/hatyan_functions.pdf doc/hatyan_functions.md'%(hatyan.__version__))

if err_html: # or err_pdf_pdoc or err_pdf_pandoc:
    raise Exception('ERROR: in documentation generation, check feedback above')
else:
    print('## succesfully finished')
