import os
import sys
import subprocess

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

##############################################################################
# Globals
##############################################################################

VERSION = '0.3'
ISRELEASED = False

##############################################################################
# Utility functions
##############################################################################


def git_version():
    """Return the git revision as a string, copied from numpy setup.py"""
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def write_version_py(filename):
    cnt = """
# THIS FILE IS GENERATED FROM TOY-ACCEL SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of numpy.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def find_packages():
    """Find all of Toy-Accel's packages.

    Credit: this code is adapted from IPython.
    https://github.com/rmcgibbo/ipython/blob/master/setupbase.py
    """
    packages = []
    for dir, subdirs, files in os.walk('src'):
        package = dir.replace(os.path.sep, '.')
        if '__init__.py' not in files:
            # not a package
            continue
        packages.append(package)
    return packages


def check_openmm_version():
    from simtk.openmm import Platform
    if not Platform.getOpenMMVersion() >= '5.1':
        raise ValueError('MSMAccelerator requires OpenMM >= 5.1')

def check_mdtraj_version():
    import pkg_resources  # part of setuptools
    err = ValueError('MSMAccelerator requires MDTraj version 0.1 or greater '
                     'https://github.com/rmcgibbo/mdtraj')
    try:
        if not pkg_resources.require("MDTraj")[0].version >= '0.1':
            raise err
    except pkg_resources.DistributionNotFound:
        raise err



##############################################################################
# Script
##############################################################################


setup_args = {
    'name': 'toy_accel',
    'version': VERSION,
    'author': 'Matthew Harrigan',
    'license': 'GPLv3',
    'url': 'https://github.com/mpharrigan/msmaccelerator-toy',
    'platforms': ['Linux', 'Mac OS X'],
    'description': 'Uses a toy system to demonstrate MSMAccelerator',
    'scripts': ['Toy.py']
}

if 'setuptools' in sys.modules:
    setup_args['zip_safe'] = False
    setup_args['install_requires'] = ['IPython>=0.12', 'pyzmq>=2.1.11',
        'pyyaml', 'sqlalchemy']


setup_args['packages'] = find_packages()
check_openmm_version()
check_mdtraj_version()
write_version_py('toy_accel/version.py')
setup(**setup_args)
