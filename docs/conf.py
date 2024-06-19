import os
import subprocess

subprocess.call('doxygen', shell=True)

extensions = ['breathe']

breathe_projects = {
    "CCEX": "../docs/xml"
}
breathe_default_project = "CCEX"

