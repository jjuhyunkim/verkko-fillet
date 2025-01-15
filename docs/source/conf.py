import sys, os
import CodeChat.CodeToRest

project = "Verkko-fillet"
copyright = "2025, NIH/NHGRI Adam Phillippy group"

version = "0.1"
release = "version 0.1"

highlight_language = "python3"
pygments_style = "sphinx"

extensions = ["CodeChat.CodeToRestSphinx", "sphinx.ext.graphviz", ]


### Options for HTML output ### 
html_theme = "alabaster"
html_static_path = CodeChat.CodeToRest.html_static_path()
html_last_updated_fmt = "%b, %d, %Y"
html_copy_source = True
html_show_sourcelink = True
html_sourcelink_suffix = ""
