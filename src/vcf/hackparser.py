import re 
import itertools

# token patterns
p_comma = r","
p_colon = r":"
p_semicolon = r";"
p_equals = r"="
p_slash = r"/"
p_bar = r"\|"
p_dot = r"\."
p_int = r"[-+]?(?:0|[1-9][0-9]*)"
p_float = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
p_str = r".*"
p_end = r'$'
p_identifier = r"[a-zA-Z][a-zA-Z0-9]*"
p_nucleotides = r"[ACGTN]+"
p_info_str = r"[^,;\n]+"
p_phase = "{slash}|{bar}".format(slash=p_slash, bar=p_bar)

def anchor(restr):
    return "^" + restr + "$"

def group(restr):
    return '(' + restr + ')'

# compiled re's
re_comma = re.compile(anchor(p_comma))
re_colon = re.compile(anchor(p_colon))
re_semicolon = re.compile(anchor(p_semicolon))
re_equals = re.compile(anchor(p_equals))
re_slash = re.compile(anchor(p_slash))
re_bar = re.compile(anchor(p_bar))
re_dot = re.compile(anchor(p_dot))
re_int = re.compile(anchor(p_int))
re_float = re.compile(anchor(p_float))
re_str = re.compile(anchor(p_str))
re_end = re.compile(anchor(p_end))
re_identifier = re.compile(anchor(p_identifier))
re_nucleotides = re.compile(anchor(p_nucleotides))
re_info_str = re.compile(anchor(p_info_str))
re_phase = re.compile(anchor(p_phase))

# python doesn't pretty print compiled regexes...
re_to_p = {
    re_comma : p_comma,
    re_colon : p_colon,
    re_semicolon : p_semicolon,
    re_equals : p_equals,
    re_slash : p_slash,
    re_bar : p_bar,
    re_dot : p_dot,
    re_int : p_int,
    re_float : p_float,
    re_str : p_str,
    re_end : p_end,
    re_identifier : p_identifier,
    re_nucleotides : p_nucleotides,
    re_info_str : p_info_str,
}

def parse(name, string):
    parser = name_to_parser[name]
    return parser(string)

def attr_restr(value_restr):
    return "(?P<attr>{attr_restr})(?:=(?P<value>{value_restr}))?".format(attr_str=r"[a-zA-Z]+", value_restr=value_restr)

typeable_as = {
        None: frozenset([int, float, bool, str]),
        int: frozenset([int, float, str]),
        float: frozenset([float, str]),
        str: frozenset([str]),
        bool: frozenset([bool]),
        }

def parse_info_attr(attr):
    attr_value = attr.split(p_equals)
    attr = attr_value[0]
    if len(attr_value) == 2:
        return (attr, parse_value(attr_value[1]))
    else:
        return (attr, True)

def parse_scalar_value(value):
    result = re_int.match(value)
    if result is not None:
        return int(value)
    result = re_float.match(value)
    if result is not None:
        return float(value)
    # return a str
    return value

def base_type(values, default=None):
    if len(values) == 0:
        return default
    types = frozenset([type(v) for v in values])
    i = iter(types)
    common_types = set(typeable_as[i.next()])
    for t in i:
        common_types.intersection_update(typeable_as[t])
    if len(common_types) == 1:
        return iter(common_types).next()
    common_types.intersection_update(types)
    if len(common_types) == 1:
        return iter(common_types).next()
    return default

def parse_value(value):
    values = value.split(p_comma)
    if len(values) != 1:
        parsed_values = [parse_scalar_value(v) for v in values]
        btype = base_type(parsed_values, default=str)
        return [btype(v) for v in parsed_values]
    return parse_scalar_value(value)

def parse_info_by_type(info):
    attrs_by_type = {
        str: [],
        int: [],
        float: [],
        bool: [],
    }
    for attr_str in info.split(p_semicolon):
        attr_pair = parse_info_attr(attr_str)  
        attrs_by_type[type(attr_pair[1])].append(attr_pair)
    return attrs_by_type

def parse_info(info):
    return dict([parse_info_attr(attr_str) for attr_str in info.split(p_semicolon)])

def ordered_alleles(ref, alts):
    if type(alts) == str:
        alts = alts.split(p_comma)
    genotypes = []
    alleles = []
    last_alleles = [ref]
    last_alleles.extend(alts)
    for x in last_alleles:
        alleles.append(x)
        for y in alleles:
            genotypes.append((y, x))
    return genotypes

def parse_dbsnp_id(string):
    return parse_none_else_string(string)

def parse_null_genotype(string):
    allele1, phase, allele2 = parse_zip_split(group(p_phase), [parse_none, parse_phase, parse_none], string)
    return ((allele1, allele2), phase)

def parse_phase(string):
    # | => True, / => False
    return match(re_phase, string, on_success=lambda r, s: True if s == "|" else False)

def parse_int(string):
    return match(re_int, string, on_success=lambda r, s: int(s))

def parse_int_list(string):
    return parse_each_split(p_comma, parse_int, string)

parse_AD = parse_int_list
parse_DP = parse_int
parse_GQ = parse_int
parse_PL = parse_int_list
def parse_GT(string):
    allele1, phase, allele2 = parse_zip_split(group(p_phase), [parse_int, parse_phase, parse_int], string)
    return ((allele1, allele2), phase)
    return tuple(parse_each_split(p_slash, parse_int, string))
def parse_nonnull_genotype(string):
    # 0/1:11,15:26:99:364,0,353
    # GT:AD:DP:GQ:PL
    GT, AD, DP, GQ, PL = re.split(p_colon, string)
    return { 'GT':parse_GT(GT), 'AD':parse_AD(AD), 'DP':parse_DP(DP), 'GQ':parse_GQ(GQ), 'PL':parse_PL(PL) }

def parse_genotype(string):
    return parse_either(parse_null_genotype, parse_nonnull_genotype, string)

def parse_genotype_format(string):
    pass

def parse_ref(string):
    return match(re_nucleotides, string)

def parse_alts(string):
    return parse_each_split(p_comma, parse_allele, string)

def parse_allele(string):
    return match(re_nucleotides, string, on_fail=lambda re, s: parse_none(s)) 

def parse_dbsnp_id(string):
    return parse_none_else_string(string)

def parse_none(string):
    return match(re_dot, string, on_success=lambda m, s: None)

def parse_none_else_string(string):
    return match(re_dot, string, on_success=lambda m, s: None, on_fail=lambda r, s: string)

def parse_each(seq, parser):
    return [parser(s) for s in seq]

def parse_zip(seq, parsers):
    return [p(s) for s, p in itertools.izip(seq, parsers)]

def parse_zip_split(regexp, parsers, string):
    return parse_zip(re.split(regexp, string), parsers)

def parse_each_split(regexp, parser, string):
    return parse_each(re.split(regexp, string), parser)

def parse_either(parser1, parser2, string):
    try:
        return parser1(string)
    except ParserException: 
        return parser2(string)

def _on_fail(regex, string):
    raise ParserException("Failed to parse {string} using {regex}".format(string=string, regex=re_to_p.get(regex, regex)))
def _on_success(match, string):
    return string
def match(regex, string, on_success=_on_success, on_fail=_on_fail):
    result = regex.match(string)
    if result is not None:
        return on_success(result, string)
    return on_fail(regex, string)

class ParserException(Exception):
    pass

name_to_parser = {
    'info' : parse_info,
    'ref' : parse_ref,
    'dbsnp_id' : parse_dbsnp_id,
    'alts' : parse_alts,
    'genotype' : parse_genotype,
    # not used
    # 'genotype_format' : parse_genotype_format,
    # 'allele' : parse_allele,
}
