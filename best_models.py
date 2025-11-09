

import os,sys,gzip,shutil,json,pickle,matplotlib,argparse
import pandas as pd
from multiprocessing import Pool

def extract_best_model(modelname,copying=True):
    ## check if the simulation is finished ##
    afdir = "../AF2_PPI/"
    rankingf = os.path.join(afdir,"%s/ranking_debug.json"%modelname)
    if os.path.isfile(rankingf):
        data = json.load(open(rankingf))
        topm = data["order"][0]
        pkl = "%s/result_%s.pkl"%(modelname,topm)
        pdb = "%s/ranked_0.pdb.gz"%modelname
        pkl = os.path.join(afdir,pkl)
        pdb = os.path.join(afdir,pdb)
        ## copying model ?
        if copying:
            shutil.copy(pdb,"%s.pdb"%modelname)
            shutil.copy(pkl,"%s.pkl"%modelname)

    ## if the prediction is not finished, check if there are any outputs ##
    else:
        max_pdb = ""
        max_pkl = ""
        max_score = 0
        for f in os.listdir(modelname):
            if f.endswith("pkl"):
                pkl = "%s/%s"%(modelname,f)
                data = pickle.load(open(pkl,"rb"))
                if "iptm" in data:
                    key = "iptm"
                elif "ptm" in data:
                    key = "ptm"
                else:
                    key = ""
                if key:
                    if data[key] > max_score:
                        max_pkl = pkl
                        max_score = data[key]
        if max_pkl == "":
            sys.stderr.write("Prediction with no outputs, quit\n")
            sys.exit(0)

        ## if there are outputs, check the first output model ##
        pkl = max_pkl
        max_name = os.path.basename(pkl)[7:-4]
        max_pdb = "%s/relaxed_%s.pdb"%(modelname,max_name)
        max_pdb = os.path.join(afdir,max_pdb)
        if os.path.isfile(max_pdb):
            pdb = max_pdb
        else:
            pdb = "%s/unrelaxed_%s.pdb"%(modelname,max_name)
            pdb = os.path.join(afdir,pdb)
        if showfigure=="Y":
            print("Not finished, plot current best model: %s"%max_name)

        ## still extract current best model ##
        shutil.copy(pdb,"%s.pdb"%modelname)
        shutil.copy(pkl,"%s.pkl"%modelname)

if __name__ == "__main__":

    df = pd.read_csv("pae_summary.csv")

    ncpu = 36
    #ncpu = 1
    if ncpu>1:
        pool = Pool(ncpu)

    for pair in df.pair:
        if ncpu==1:
            extract_best_model(pair)
        else:
            pool.apply_async(extract_best_model,args=(pair,))
    if ncpu>1:
        pool.close()
        pool.join()