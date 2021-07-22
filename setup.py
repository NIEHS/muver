#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0,<=7.0',
    'numpy>=1.12.0,<=1.16.6',
    'scipy>=0.18.1,<=1.2.3',
    'matplotlib>=2.0.0,<=2.2.3',
    'regex>=2017.6.23,<=2019.12.20',
]

setup_requirements = [
    'pytest-runner<=5.2',
]

test_requirements = [
    'pytest<=4.6.11',
]

setup(
    name='muver',
    version='1.2.4',
    description="SNP and indel caller for mutation accumulation experiments",
    long_description=readme + '\n\n' + history,
    author="Christopher Andrew Lavender, Adam Burkholder",
    author_email='adam.burkholder@nih.gov',
    url='https://github.com/NIEHS/muver',
    packages=find_packages(include=['muver', 'muver.wrappers']),
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
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
    data_files=[('config', ['paths.cfg'])],
)
