import os
import sys

sys.path.insert(0, os.path.abspath('.'))

extensions = ['breathe']

breathe_projects = {
    "CCEX": "docs/xml"
}
breathe_default_project = "CCEX"

