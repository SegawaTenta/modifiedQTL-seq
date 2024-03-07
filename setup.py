#!/usr/bin/env python3

from setuptools import setup, find_packages
from modified.__init__ import __version__

setup(
    name="modified QTL-seq",
    version='{}'.format(__version__),
    description='modified QTL-seq: pipeline to identify causative genomic region responsible for a phenotype',
    url='https://github.com/SegawaTenta/modifiedQTL-seq',
    author='Tenta Segawa',
    packages=['mqtlseq'],
    entry_points={
        'console_scripts': [
            'mqtlseq=mqtlseq.mqtlseq:main'
            'mqtlseq_mpileup=mqtlseq.mqtlseq_mpileup:main'
            'mqtlseq_snpindex=mqtlseq.mqtlseq_snpindex:main'
            'mqtlseq_plot=mqtlseq.mqtlseq_plot:main'
        ]
    }
)
