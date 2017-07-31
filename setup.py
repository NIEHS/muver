#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    'numpy>=1.12.0',
    'scipy>=0.18.1',
    'matplotlib>=2.0.0',
    'regex>=2017.6.23',
    # TODO: put package requirements here
]

setup_requirements = [
    'pytest-runner',
    # TODO(lavenderca): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'pytest',
    # TODO: put package test requirements here
]

setup(
    name='muver',
    version='0.1.0',
    description="SNP and indel caller for mutation accumulation experiments",
    long_description=readme + '\n\n' + history,
    author="Christopher Andrew Lavender",
    author_email='christopher.lavender@nih.gov',
    url='https://github.com/lavenderca/muver',
    packages=find_packages(include=['muver']),
    entry_points={
        'console_scripts': [
            'muver=muver.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='muver',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
    data_files=[('config', ['paths.cfg'])],
)
