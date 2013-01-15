Installation:
The python vcf module uses the yapps library (http://pypi.python.org/pypi/Yapps2/2.0.4) to generate a parser 
(vcfparser.py) for parsing vcf fields.

To generate the parser, simply run:

make -f vcf.mk

Then, using the parser in your code is as easy as:

import vcf

vcf.parse('genotype', '0/0:23,0,0:23:60:0,60,709,60,709,709')
>>> 
    {'AD': [23, 0, 0],
     'DP': 23,
     'GQ': 60,
     'GT': ((0, 0), False),
     'PL': [0, 60, 709, 60, 709, 709]}
vcf.parse('info', 'AC=1;AF=0.042;AN=24;BaseQRankSum=1.008;DP=280;Dels=0.00;FS=2.163;HaplotypeScore=1.0714;InbreedingCoeff=-0.0452;MLEAC=1;MLEAF=0.042;MQ=59.92;MQ0=0;MQRankSum=-0.875;QD=8.58;ReadPosRankSum=0.314;SB=-1.641e-02;VQSLOD=5.1338;culprit=QD')
>>>
    {'AC': 1,
     'AF': 0.042,
     'AN': 24,
     'BaseQRankSum': 1.008,
     'DP': 280,
     'Dels': 0.0,
     'FS': 2.163,
     'HaplotypeScore': 1.0714,
     'InbreedingCoeff': -0.0452,
     'MLEAC': 1,
     'MLEAF': 0.042,
     'MQ': 59.92,
     'MQ0': 0,
     'MQRankSum': -0.875,
     'QD': 8.58,
     'ReadPosRankSum': 0.314,
     'SB': -0.01641,
     'VQSLOD': 5.1338,
     'culprit': 'QD'}

To see other things you can parse, look at the rules in vcfparser.g that look like:
rule info: (
    ...
) END
