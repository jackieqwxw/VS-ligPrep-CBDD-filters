#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt, MolLogP,  NumHDonors, NumHAcceptors
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
import multiprocessing as mp
from multiprocessing import Pool
import time
import pandas as pd
import os
import json
from docopt import docopt
import pkg_resources
from rdkit.Chem.Lipinski import RingCount

cmd_str = """Usage:
new_filters.py filter --in INPUT_FILE --prefix PREFIX [--rules RULES_FILE_NAME] [--alerts ALERT_FILE_NAME][--np NUM_CORES]
new_filters.py template --out TEMPLATE_FILE [--rules RULES_FILE_NAME]

Options:
--in INPUT_FILE input file name
--prefix PREFIX prefix for output file names
--rules RULES_FILE_NAME name of the rules JSON file
--alerts ALERTS_FILE_NAME name of the structural alerts file
--np NUM_CORES the number of cpu cores to use (default is all)
--out TEMPLATE_FILE parameter template file name
"""

## Number of cycles > 0
def GetRingSystems(mol, includeSpiro=True):         
    """
    list of set of ring systems in each mol
    """
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems

## Number of chiral centers < 3
def GetNumChiralCenter(mol):
    chi_cen_len=len(Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=False,useLegacyImplementation=True))
    return chi_cen_len

## Number of carboxyl groups <= 1
def GetNumCarboxylGroups(mol):
    return Chem.Fragments.fr_COO2(mol)

## Number of amino N with aliphatic neighbors <=2
def GetNumAminoNitrogen(mol): 
    atoms=[x for x in mol.GetAtoms()] 
    ind=[x.GetIdx() for x in atoms]
    atom_num=[x.GetAtomicNum() for x in atoms]
    atom_hyb=[x.GetHybridization() for x in atoms]
    comb=list(zip(atom_num, atom_hyb))
    Natoms=[(i[1]) for _,i in enumerate(comb) if i[0]==7]
    count=len([x for x in Natoms if x == Chem.rdchem.HybridizationType.SP3])
    comb=list(zip(ind,atom_num))
    a=[i for i,j in comb if j==7]
    b=[i for i,j in comb if j==16]
    c=[[i+1, i-1, i+2, i-2, i+3, i-3, i+4, i-4, i+5, i-5] for i,j in comb if j==7]
    m=[]
    n=[]
    for i in b:
        for j in c:
            for k in j:
                if i==k:
                    m.append(i)
                else:
                    None
    count=max(count-len(m),0)
    return count

## Number of negative ionizable + Number of positive ionizable <= 1
def GetNumCarboxylAndAmino(mol):
    car_count = GetNumCarboxylGroups(mol)
    aminoN_count = GetNumAminoNitrogen(mol)
    tot = car_count + aminoN_count
    return tot

## Number of fused rings in a single system < 4 
def GetNumFusedRing(mol):
    ri = mol.GetRingInfo()
    d=[]
    for ring in ri.AtomRings():
        d.append(set(ring))
    total_rings=len(d)
    fused_bonds_ring=0
    fused_atoms=[]
    fu_rings=0
    count_fused_bonds=0
    from itertools import combinations
    for items in list(combinations(d, 2)):
        if items[0].intersection(items[1]):
            fused_atoms.append(items[0].union(items[1]))
            count_fused_bonds+=1
            if len(fused_atoms)>1:
                for items in list(combinations(fused_atoms, 2)):
                    if items[0].intersection(items[1]):
                        fused_bonds_ring=+1
    return fused_bonds_ring  

## Number of acyclic N-N (or N=N) bonds = 0 (hydrazone)
def GetNumHydrazone(mol):
    i = []
    for bond in mol.GetBonds():
        btype = bond.GetBondType()
        b1AtNum = bond.GetBeginAtom().GetAtomicNum()
        b2AtNum = bond.GetEndAtom().GetAtomicNum()
        b1AtRing = bond.GetBeginAtom().IsInRing()
        b2AtRing = bond.GetEndAtom().IsInRing()
        if b1AtNum == 7 and b2AtNum == 7 and b1AtRing == False and b2AtRing == False or btype == 'DOUBLE' or btype == 'SINGLE' :
            i.append(btype)
    i_len = len(i)
    return i_len

## Number of long aliphatic chains (Nbonds>=6) = 0
def aliphatic_atoms(mol): 
    """Gives the generator of list of aliphatic chains if present in each compounds"""
    rot_atom_pairs = list(mol.GetSubstructMatches(Chem.MolFromSmarts("[R0;D2]")))
    l=[list(x) for x in rot_atom_pairs]
    f=[(i[0]) for i in l]
    import itertools
    for i, j in itertools.groupby(enumerate(f), lambda x: x[1] - x[0]):
        j = list(j)
        start = j[0][1]
        length = len(j)
        if length == 1:
            yield (start-start)
        else:
            yield ((start+length)-start)
def connect_aa(mol):
    ali_atoms=[]
    ali_atoms.append(list(aliphatic_atoms(mol)))        
    return ali_atoms
def Aliphatic_c(ali_atoms):
    ali_c=[]
    for i in ali_atoms:
        if i and max(i):
            ali_c.append(max(i))
        else:
            ali_c.append(0)
    return ali_c
def aliphatic_atom_count(mol):
    """overall program including aliphatic_atoms(),connect_aa(),Aliphatic_c(),aliphatic_atom_count()
    to count the longest aliphatic chain in a given molecule"""
    result=connect_aa(mol)
    ali_c= Aliphatic_c(list(result))
    ali_c_ele = ali_c[0]
    return ali_c_ele


def read_rules(rules_file_name):
    """
    Read rules from a JSON file
    :param rules_file_name: JSON file name
    :return: dictionary corresponding to the contents of the JSON file
    """
    with open(rules_file_name) as json_file:
        try:
            rules_dict = json.load(json_file)
            return rules_dict
        except json.JSONDecodeError:
            print(f"Error parsing JSON file {rules_file_name}")
            sys.exit(1)


def write_rules(rule_dict, file_name):
    """
    Write configuration to a JSON file
    :param rule_dict: dictionary with rules
    :param file_name: JSON file name
    :return: None
    """
    ofs = open(file_name, "w")
    ofs.write(json.dumps(rule_dict, indent=4, sort_keys=True))
    print(f"Wrote rules to {file_name}")
    ofs.close()


def default_rule_template(alert_list, file_name):
    """
    Build a default rules template
    :param alert_list: list of alert set names
    :param file_name: output file name
    :return: None
    """
    default_rule_dict = {
        "MW": [300, 600],
        "HBD": [0, 5],
        "HBA": [0, 12],
        "LogP": [1.00001, 5.99999],
        "RotBond": [0, 12],
        "Hydrazone": [0, 0],
        "FusedRing": [0, 3],
        "Carboxyl": [0, 1],
        "Amino": [0, 2],
        "CarboxAmino": [0, 1],
        "Chiral": [0, 2],
        "Cycle": [1, 9999],
        "AliphChain": [0, 5]
    }
    for rule_name in alert_list:
        if rule_name == "Inpharmatica":
            default_rule_dict["Rule_" + rule_name] = True
        else:
            default_rule_dict["Rule_" + rule_name] = False
    write_rules(default_rule_dict, file_name)


def get_config_file(file_name, environment_variable):
    """
    Read a configuration file, first look for the file, if you can't find
    it there, look in the directory pointed to by environment_variable
    :param file_name: the configuration file
    :param environment_variable: the environment variable
    :return: the file name or file_path if it exists otherwise exit
    """
    if os.path.exists(file_name):
        return file_name
    else:
        config_dir = os.environ.get(environment_variable)
        if config_dir:
            config_file_path = os.path.join(os.path.sep, config_dir, file_name)
            if os.path.exists(config_file_path):
                return config_file_path

    error_list = [f"Could not file {file_name}"]
    if config_dir:
        err_str = f"Could not find {config_file_path} based on the {environment_variable}" + \
                  "environment variable"
        error_list.append(err_str)
    error_list.append(f"Please check {file_name} exists")
    error_list.append(f"Or in the directory pointed to by the {environment_variable} environment variable")
    print("\n".join(error_list))
    sys.exit(1)


class RDFilters:
    def __init__(self, rules_file_name):
        good_name = get_config_file(rules_file_name, "FILTER_RULES_DIR")
        self.rule_df = pd.read_csv(good_name)
        # make sure there wasn't a blank line introduced
        self.rule_df = self.rule_df.dropna()
        self.rule_list = []

    def build_rule_list(self, alert_name_list):
        """
        Read the alerts csv file and select the rule sets defined in alert_name_list
        :param alert_name_list: list of alert sets to use
        :return:
        """
        self.rule_df = self.rule_df[self.rule_df.rule_set_name.isin(alert_name_list)]
        tmp_rule_list = self.rule_df[["rule_id", "smarts", "max", "description"]].values.tolist()
        for rule_id, smarts, max_val, desc in tmp_rule_list:
            smarts_mol = Chem.MolFromSmarts(smarts)
            if smarts_mol:
                self.rule_list.append([smarts_mol, max_val, desc])
            else:
                print(f"Error parsing SMARTS for rule {rule_id}", file=sys.stderr)

    def get_alert_sets(self):
        """
        :return: a list of unique rule set names
        """
        return self.rule_df.rule_set_name.unique()

    def evaluate(self, lst_in):
        """
        Evaluate structure alerts on a list of SMILES
        :param lst_in: input list of [SMILES, Name]
        :return: list of alerts matched or "OK"
        """
        smiles, name = lst_in
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [smiles, name, 'INVALID', -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999]
        desc_list = [MolWt(mol), NumHDonors(mol), NumHAcceptors(mol), MolLogP(mol), CalcNumRotatableBonds(mol),  GetNumHydrazone(mol), GetNumFusedRing(mol), GetNumCarboxylGroups(mol), GetNumAminoNitrogen(mol), GetNumCarboxylAndAmino(mol), GetNumChiralCenter(mol), RingCount(mol), aliphatic_atom_count(mol)]
        for row in self.rule_list:
            patt, max_val, desc = row
            if len(mol.GetSubstructMatches(patt)) > max_val:
                return [smiles, name] + [desc + " > %d" % (max_val)] + desc_list
        return [smiles, name] + ["OK"] + desc_list


def main():
    cmd_input = docopt(cmd_str)
    alert_file_name = cmd_input.get("--alerts") or pkg_resources.resource_filename('new_filters.py',
                                                                                   "data/alert_collection.csv")
    rf = RDFilters(alert_file_name)

    if cmd_input.get("template"):
        template_output_file = cmd_input.get("--out")
        default_rule_template(rf.get_alert_sets(), template_output_file)

    elif cmd_input.get("filter"):
        input_file_name = cmd_input.get("--in")
        rules_file_name = cmd_input.get("--rules") or pkg_resources.resource_filename('new_filters.py', "data/rules.json")
        rules_file_path = get_config_file(rules_file_name, "FILTER_RULES_DATA")
        prefix_name = cmd_input.get("--prefix")
        num_cores = cmd_input.get("--np") or mp.cpu_count()
        num_cores = int(num_cores)

        print("using %d cores" % num_cores, file=sys.stderr)
        start_time = time.time()
        p = Pool(num_cores)
        input_data = [x.split() for x in open(input_file_name)]
        input_data = [x for x in input_data if len(x) == 2]
        rule_dict = read_rules(rules_file_path)

        rule_list = [x.replace("Rule_", "") for x in rule_dict.keys() if x.startswith("Rule") and rule_dict[x]]
        rule_str = " and ".join(rule_list)
        print(f"Using alerts from {rule_str}", file=sys.stderr)
        rf.build_rule_list(rule_list)
        res = list(p.map(rf.evaluate, input_data))
        df = pd.DataFrame(res, columns=["SMILES", "NAME", "FILTER", "MW", "HBD", "HBA", "LogP", "RotBond","Hydrazone", "FusedRing", "Carboxyl", "Amino", "CarboxAmino", "Chiral", "Cycle", "AliphChain"])
        df_ok = df[
            (df.FILTER == "OK") &
            df.MW.between(*rule_dict["MW"]) &
            df.HBD.between(*rule_dict["HBD"]) &
            df.HBA.between(*rule_dict["HBA"]) & 
            df.LogP.between(*rule_dict["LogP"]) & 
            df.RotBond.between(*rule_dict["RotBond"]) & 
            df.Hydrazone.between(*rule_dict["Hydrazone"]) &
            df.FusedRing.between(*rule_dict["FusedRing"]) &
            df.Carboxyl.between(*rule_dict["Carboxyl"]) &
            df.Amino.between(*rule_dict["Amino"]) &
            df.CarboxAmino.between(*rule_dict["CarboxAmino"]) &
            df.Chiral.between(*rule_dict["Chiral"]) & 
            df.Cycle.between(*rule_dict["Cycle"]) &
            df.AliphChain.between(*rule_dict["AliphChain"]) 
            ]
        output_smiles_file = prefix_name + ".smi"
        output_csv_file = prefix_name + ".csv"
        df_ok[["SMILES", "NAME"]].to_csv(f"{output_smiles_file}", sep=" ", index=False, header=False)
        print(f"Wrote SMILES for molecules passing filters to {output_smiles_file}", file=sys.stderr)
        df.to_csv(f"{prefix_name}.csv", index=False)
        print(f"Wrote detailed data to {output_csv_file}", file=sys.stderr)

        num_input_rows = df.shape[0]
        num_output_rows = df_ok.shape[0]
        fraction_passed = "%.1f" % (num_output_rows / num_input_rows * 100.0)
        print(f"{num_output_rows} of {num_input_rows} passed filters {fraction_passed}%", file=sys.stderr)
        elapsed_time = "%.2f" % (time.time() - start_time)
        print(f"Elapsed time {elapsed_time} seconds", file=sys.stderr)


if __name__ == "__main__":
    main()
