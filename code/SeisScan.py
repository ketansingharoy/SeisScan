import os
import sys

#--------- Generate important path ----------------
code_path = os.path.dirname(os.path.abspath(__file__))

# print(code_path)

lib_path = os.path.join(code_path, 'lib')
script_path = os.path.join(code_path, 'script')
utils_path = os.path.join(code_path, 'utils')

# lib_path = os.path.join(code_path, '..', 'lib')
# lib_path = os.path.realpath(lib_path)

# script_path = os.path.join(code_path, '..', 'script')
# script_path = os.path.realpath(script_path)

# utils_path = os.path.join(code_path, '..', 'utils')
# utils_path = os.path.realpath(utils_path)


#---------- Set as system path ---------------
sys.path.insert(1, lib_path)
sys.path.insert(2, script_path)
sys.path.insert(3, utils_path)

#---------- import function from modules in utila ---------------
from utils import db, prs

make_sds_from_fdsn = db.make_sds_from_fdsn
read_sds = db.read_sds
get_inv_from_db = db.get_inv_from_db
prs = prs.prs

#----------- import functions from modules in script --------------
from script import peak_cross_correlation, local_similarity, backprojection

pcc = peak_cross_correlation.do_pcc
do_ls = local_similarity.do_ls

prepare_traveltime_lookup_table = backprojection.prepare_traveltime_lookup_table
do_bp = backprojection.do_bp

#----------- remove from system path
# sys.path.remove(lib_path)
# sys.path.remove(script_path)
# sys.path.remove(utils_path)
