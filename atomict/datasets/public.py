# This is a public s3 bucket

ROOT_URL = 'https://atomict-raw-datasets.s3.us-west-2.amazonaws.com/'

MAPPINGS = {
    'bace': f'{ROOT_URL}bace.csv.gz',
    'BBBP': f'{ROOT_URL}BBBP.csv.gz',
    'clintox': f'{ROOT_URL}clintox.csv.gz',
    'delaney': f'{ROOT_URL}delaney-processed.csv.gz',
    'freesolv': f'{ROOT_URL}freesolv.csv.gz',
    'gdb7': f'{ROOT_URL}gdb7.tar.gz',
    'gdb8': f'{ROOT_URL}gdb8.tar.gz',
    'gdb9': f'{ROOT_URL}gdb9.tar.gz',
    'hiv': f'{ROOT_URL}HIV.csv.gz',
    'hopv': f'{ROOT_URL}hopv.csv.gz',
    'hppb': f'{ROOT_URL}hppb.csv.gz',
    'lipophilicity': f'{ROOT_URL}lipophilicity.csv.gz',
    'muv': f'{ROOT_URL}muv.csv.gz',
    'nci': f'{ROOT_URL}nci_unique.csv.gz',
    'pcba': f'{ROOT_URL}pcba.csv.gz',
    'pdbbind_refined': f'{ROOT_URL}pdbbind_v2019_refined.tar.gz',
    'pdbbind': f'{ROOT_URL}PDBBind.zip',
    'pfas_4k': f'{ROOT_URL}pfas_4k.csv.gz',
    'pfas_priority': f'{ROOT_URL}pfas_priority.csv.gz',
    'RGD1': f'{ROOT_URL}RGD1_CHNO.h5.gz',
    'sider': f'{ROOT_URL}sider.csv.gz',
    'thermosol': f'{ROOT_URL}thermosol.csv.gz',
    'tox21': f'{ROOT_URL}tox21.csv.gz',
    'toxcast': f'{ROOT_URL}toxcast_data.csv.gz',
    'USPTO': f'{ROOT_URL}USPTO-FULL.csv.gz',
    'zinc15': f'{ROOT_URL}zinc15_270M_2D.csv.gz',
}