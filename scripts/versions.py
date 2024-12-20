#!/usr/bin/env python3
"""
This script uses pip to print out versions of all installed packages.
"""

import pip

installed_packages = pip.get_installed_distributions()
installed_packages_list = sorted(
    ["%s==%s" % (i.key, i.version) for i in installed_packages]
)

print(installed_packages_list)
