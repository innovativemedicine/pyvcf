import vcf
import csv
import fileinput
import sys

# for vc_group in vcf_file("data/Project_Charames_TruSight_Cancer_Panel.130117.raw.snps.genome_summary.csv"):
#     vc_group.columns['chromosome']
#     vc_group.vc_group_allele # list
#     for vc in vc_group.vc:
#         vc.columns['genotype_quality']
#         vc.vc_allele # list
#         vc.vc_allele[0].columns
#         vc.vc_genotype # list
def vcf_file(filepath=None, input=None, delim=",", quote='"'):
    if filepath is not None:
        input = open(filepath, 'rb')
    csv_input = csv.reader(input, delimiter=delim, quotechar=quote)

    # for each vc_group
    for row in csv_input:

        row = [None if f == '' else f for f in row] 

        info = vcf.parse('info', row[35 - 1])
        vc_group_columns = {
            'chromosome'  : row[22 - 1],
            'start_posn'  : row[23 - 1],
            'end_posn'    : row[24 - 1],
            'ref'         : vcf.parse('ref', row[25 - 1]),
            'dbsnp_id'    : vcf.parse('dbsnp_id', row[30 - 1]),

            # 'genotype_format' : row[36 - 1],
            'quality'         : row[33 - 1],
            'filter'          : row[34 - 1],

            # annovar columns
            'otherinfo'               : row[27 - 1],
            'func'                    : row[1 - 1],
            'gene'                    : row[2 - 1],
            'exonicfunc'              : row[3 - 1],
            'aachange'                : row[4 - 1],
            'conserved'               : row[5 - 1],
            '1000g2011may_all'        : row[8 - 1],
            'dbsnp135'                : row[9 - 1],
            'ljb_phylop_pred'         : row[12 - 1],
            'ljb_sift_pred'           : row[14 - 1],
            'ljb_polyphen2_pred'      : row[16 - 1],
            'ljb_lrt_pred'            : row[18 - 1],
            'ljb_mutationtaster_pred' : row[20 - 1],

            'ljb_gerppp'              : row[21 - 1],
            'segdup'                  : row[6 - 1],
            'esp5400_all'             : row[7 - 1],
            'avsift'                  : row[10 - 1],
            'ljb_phylop'              : row[11 - 1],
            'ljb_sift'                : row[13 - 1],
            'ljb_polyphen2'           : row[15 - 1],
            'ljb_lrt'                 : row[17 - 1],
            'ljb_mutationtaster'      : row[19 - 1],

            # vc_group_info columns
            # 'info_source'       : row[36 - 1],
            'ds'                : info.get('DS', False),
            'inbreeding_coeff'  : info.get('InbreedingCoeff'),
            'base_q_rank_sum'   : info.get('BaseQRankSum'),
            'mq_rank_sum'       : info.get('MQRankSum'),
            'read_pos_rank_sum' : info.get('ReadPosRankSum'),
            'dels'              : info.get('Dels'),
            'fs'                : info.get('FS'),
            'haplotype_score'   : info.get('HaplotypeScore'),
            'mq'                : info.get('MQ'),
            'qd'                : info.get('QD'),
            'sb'                : info.get('SB'),
            'vqslod'            : info.get('VQSLOD'),
            'an'                : info.get('AN'),
            'dp'                : info.get('DP'),
            'mq0'               : info.get('MQ0'),
            'culprit'           : info.get('culprit'),
        }
        vc_group_table = vc_group(vc_group_columns)

        alts = vcf.parse('alts', row[32 - 1])

        # for each vc_group_allele in (vc x alt alleles in vc_group)
        vc_group_allele_fields = [
            alts,
            # vc_group_allele_info
            get_list(info, 'AF'),
            get_list(info, 'MLEAF'),
            get_list(info, 'AC'),
            get_list(info, 'MLEAC'),
        ]

        vc_group_allele_columns = [{
            # 'vc_group_id' : vc_group_id,
            'allele'        : allele,
            'af'            : af,
            'mle_af'        : mle_af,
            'ac'            : ac,
            'mle_ac'        : mle_ac,
        } for allele, af, mle_af, ac, mle_ac in arity_zip(vc_group_allele_fields, table='vc_group', key="alt alleles in vc_group")]
        add_columns(vc_group_allele_columns, vc_group_table.vc_group_allele)

        ref_and_alts = as_list(vc_group_columns['ref']) + alts

        # for each vc in vc_group
        for genotype in [vcf.parse('genotype', row[gf]) for gf in xrange(37 - 1, 48)]:
            # vc_columns['genotype_source'] = row[gf]
            vc_columns = {
                # 'vc_group_id' : vc_group_table.lastrowid,
                'zygosity'    : row[27 - 1],
            }

            patient_columns = {
            }

            # vc_columns['patient_id'] = patient_table.lastrowid
            if not (type(genotype) == tuple and genotype[0] == (None, None)):
                ((allele1_idx, allele2_idx), vc_columns['phased']) = genotype['GT'] 
                vc_columns['allele1'] = ref_and_alts[allele1_idx]
                vc_columns['allele2'] = ref_and_alts[allele2_idx]
                vc_columns['read_depth'] = genotype.get('DP')
                vc_columns['genotype_quality'] = genotype.get('GQ')
                vc_table = vc(vc_columns)
                vc_group_table.vc.append(vc_table)
                
                # for each vc_genotype in (alleles in vc_group x alleles in vc_group x vc)
                vc_genotype_fields = [
                    vcf.ordered_alleles(vc_group_columns['ref'], alts), 
                    as_list(genotype.get('PL')),
                ]
                vc_genotype_columns = [{
                    # 'vc_id': vc_id,
                    'allele1': vc_genotype_allele1,
                    'allele2': vc_genotype_allele2,
                    'phred_likelihood': phred_likelihood,
                } for (vc_genotype_allele1, vc_genotype_allele2), phred_likelihood in arity_zip(vc_genotype_fields, table='vc_genotype', key="biallelic genotypes in vc_group")]
                add_columns(vc_genotype_columns, vc_table.vc_genotype)

                # for each vc_allele in (vc x alleles in ref, alts)
                vc_allele_fields = [
                    ref_and_alts,
                    get_list(genotype, 'AD'),
                ]
                vc_allele_columns = [{
                    # 'vc_id': vc_id,
                    'allele': allele,
                    'allelic_depth': allelic_depth,
                } for allele, allelic_depth in arity_zip(vc_allele_fields, table='vc_allele', key="ref and alt alleles in vc_group")]
                add_columns(vc_allele_columns, vc_table.vc_allele)

        yield vc_group_table
    input.close() 

def arity_zip(args, error=None, table=None, key=None):
    if error is None:
        error = "Number of {table} columns don't all match the number of {key}; " + \
                "skipping insertion into {table} at line {lineno}"
    return check_arity_zip(args, error % { 'lineno':1, 'table':table if type(table) == str else table.name, 'key':key })

def as_list(x):
    return [x] if type(x) != list else x

def get_list(dic, attr):
    return as_list(dic.get(attr, []))

def check_arity_zip(args, error=None):
    if not all(len(f) == len(args[0]) for f in args):
        if error is not None:
            print >> sys.stderr, error
        return []
    return zip(*args)


def add_columns(columns, table_list):
    if type(columns) == list:
        for c in columns:
            table_list.append(table(c))
    else:
        table_list.append(table(columns))

class table(object):
    def __init__(self, columns=None):
        self.columns = columns

class vc_group(table):
    def __init__(self, columns=None):
        super(vc_group, self).__init__(columns)
        self.vc = []
        self.vc_group_allele = []

class vc(table):
    def __init__(self, columns=None):
        super(vc, self).__init__(columns)
        self.vc_allele = []
        self.vc_genotype = []
