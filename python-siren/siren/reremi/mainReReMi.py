#!/usr/bin/python

import sys, re, os.path, datetime, glob
import numpy
import tempfile

import codecs

from toolLog import Log
from classPackage import Package, saveAsPackage, writeRedescriptions, getPrintParams
from classData import Data, NA_str_def
from classRedescription import Redescription
from classBatch import Batch
from classConstraints import Constraints
from classPreferencesManager import PreferencesReader, getPM
from classMiner import instMiner, StatsMiner
from classQuery import Op
## from classMinerFolds import instMiner, StatsMiner
from classQuery import Query

from classRndFactory import RndFactory

import pdb

## from codeRRM import RedModel

delims_dict = {"(auto)": None,
               "TAB": '\t',
               "SPC": ' '}

def loadAll(arguments=[], conf_defs=None):
    pm = getPM(conf_defs)

    exec_folder = os.path.dirname(os.path.abspath(__file__))
    src_folder = exec_folder
    
    pack_filename = None
    config_filename = None
    tmp_dir = None
    params = None
    reds = None
    options_args = arguments[1:]

    if len(arguments) > 1:
        if arguments[1] == "--config":
            print PreferencesReader(pm).dispParameters(None, True, True, True)
            sys.exit(2)
        if os.path.isfile(arguments[1]):
            if os.path.splitext(arguments[1])[1] == Package.DEFAULT_EXT:
                pack_filename = arguments[1]
                if len(arguments) > 2 and os.path.isfile(arguments[2]):
                    config_filename = arguments[2]
                    options_args = arguments[3:]
                else:
                    options_args = arguments[2:]
            else:
                config_filename = arguments[1]
                options_args = arguments[2:]

    if pack_filename is not None:
        src_folder = os.path.dirname(os.path.abspath(pack_filename))

        package = Package(pack_filename)
        elements_read = package.read(pm)        
        data = elements_read.get("data", None)
        reds = elements_read.get("reds", None)
        params = elements_read.get("preferences", None)
        tmp_dir = package.getTmpDir()
        
    elif config_filename is not None:
        src_folder = os.path.dirname(os.path.abspath(config_filename))

    queries_second = None
    try:
        params = PreferencesReader(pm).getParameters(config_filename, options_args, params)        
    except AttributeError:
        queries_second = config_filename
        
    if params is None:
        print 'ReReMi redescription mining\nusage: "%s [package] [config_file]"' % arguments[0]
        print '(Type "%s --config" to generate a default configuration file' % arguments[0]
        sys.exit(2)
    
    params_l = turnToDict(params)
    filenames = prepareFilenames(params_l, tmp_dir, src_folder)
    if queries_second is not None:
        filenames["queries_second"] = queries_second
    logger = Log(params_l['verbosity'], filenames["logfile"])

    if pack_filename is None:
        data = Data([filenames["LHS_data"], filenames["RHS_data"]]+filenames["add_info"], filenames["style_data"])
    logger.printL(2, data, "log")

    if pack_filename is not None:
        filenames["package"] = os.path.abspath(pack_filename)
    print filenames
    return {"params": params, "data": data, "logger": logger,
            "filenames": filenames, "reds": reds, "pm": pm}

def turnToDict(params):
    params_l = {}
    for k, v in  params.items():
        params_l[k] = v["data"]
    return params_l

def getDataAddInfo(params_l={}):
    return [{}, params_l.get('NA_str', NA_str_def)]

def prepareFilenames(params_l, tmp_dir=None, src_folder=None):
    filenames = {"queries": "-",
                 "style_data": "csv",
                 "add_info": getDataAddInfo(params_l)
                 }

    if 'delim_in' in params_l:
        dl = delims_dict.get(params_l['delim_in'], params_l['delim_in'])
        if dl is not None:
            filenames["add_info"][0]["delimiter"] = dl
    
    for p in ['result_rep', 'data_rep']:
        if p not in params_l:
            params_l[p] = ""
        if src_folder is not None and re.match("./", params_l[p]):
            params_l[p] = src_folder+params_l[p][1:]
        elif params_l[p] == "__TMP_DIR__":
            if tmp_dir is None:
                tmp_dir = tempfile.mkdtemp(prefix='ReReMi')
            params_l[p] = tmp_dir + "/"
        elif sys.platform != 'darwin':
            params_l[p] = re.sub("~", os.path.expanduser("~"), params_l[p])

    ### Make data file names
    filenames["LHS_data"] = ""
    if len(params_l.get("LHS_data", "")) != 0:
        filenames["LHS_data"] = params_l['LHS_data']
    elif len(params_l.get('data_l', "")) != 0:
        filenames["LHS_data"] = params_l['data_rep']+params_l['data_l']+params_l.get('ext_l', "")

    filenames["RHS_data"] = ""
    if len(params_l.get("RHS_data", "")) != 0 :
        filenames["RHS_data"] = params_l['RHS_data']
    elif len(params_l.get('data_r', "")) != 0:
        filenames["RHS_data"] = params_l['data_rep']+params_l['data_r']+params_l.get('ext_r', "")

    if len(params_l.get("trait_data", "")) != 0 :
        filenames["traits_data"] = params_l['traits_data']
    elif len(params_l.get('data_t', "")) != 0:
        filenames["traits_data"] = params_l['data_rep']+params_l['data_t']+params_l.get('ext_t', "")

        
    if os.path.splitext(filenames["LHS_data"])[1] != ".csv" or os.path.splitext(filenames["RHS_data"])[1] != ".csv":
        filenames["style_data"] = "multiple"
        filenames["add_info"] = []

    ### Make queries file names
    if len(params_l.get("queries_file", "")) != 0 :
        filenames["queries"] = params_l["queries_file"]
    elif params_l.get('out_base', "-") != "-"  and len(params_l['out_base']) > 0 and len(params_l.get('ext_queries', ".queries")) > 0:
        filenames["queries"] = params_l['result_rep']+params_l['out_base']+params_l.get('ext_queries', ".queries")

    if filenames["queries"] != "-":
        if not os.path.isfile(filenames["queries"]):
            try:
                tfs = open(filenames["queries"], "a")
                tfs.close()
            except IOError:
                print "Queries output file not writable, using stdout instead..."
                filenames["queries"] = "-"
    parts = filenames["queries"].split(".")
    basis = ".".join(parts[:-1])
    filenames["basis"] = basis

    ### Make named queries file name
    if filenames["queries"] != "-" and params_l.get("queries_named_file", "") == "+":
        filenames["queries_named"] = basis+"_named."+parts[-1]
    elif len(params_l.get("queries_named_file", "")) > 0:
        filenames["queries_named"] = params_l["queries_named_file"]

    ### Make support file name
    if filenames["queries"] != "-" and params_l.get("support_file", "") == "+" and len(params_l.get('ext_support', "")) > 0:
        filenames["support"] = basis+params_l['ext_support']
    elif len(params_l.get("support_file", "")) > 0:
        filenames["support"] = params_l["support_file"]

    ### Make log file name
    if filenames["queries"] != "-" and params_l.get('logfile', "") == "+" and len(params_l.get('ext_log', "")) > 0:
        filenames["logfile"] = basis+params_l['ext_log']
    elif len(params_l.get('logfile', "")) > 0:
        filenames["logfile"] = params_l['logfile']

    if len(params_l.get("series_id", "")) > 0:
        for k in filenames.keys():
            if type(filenames[k]) is str:
                filenames[k] = filenames[k].replace("__SID__", params_l["series_id"])

    return filenames

def outputResults(filenames, results_batch, data=None, with_headers=True, mode="w", data_recompute=None):
    rp = Redescription.getRP()
    modifiers, modifiers_recompute = {}, {}
    if data is not None:
        modifiers = rp.getModifiersForData(data)
    if data_recompute is not None:
        modifiers_recompute = rp.getModifiersForData(data_recompute)
    fstyle = "basic"
    
    header_recompute = ""
    if data_recompute is not None:
        fields_recompute = rp.getListFields("stats", modifiers_recompute)
        header_recompute = rp.dispHeaderFields(fields_recompute) + "\tacc_diff"

    filesfp = {"queries": None, "queries_named": None, "support": None}
    if filenames["queries"] == "-":
        filesfp["queries"] = sys.stdout
    else:
        filesfp["queries"] = open(filenames["queries"], mode)
    all_fields = rp.getListFields(fstyle, modifiers)
    if with_headers:
        filesfp["queries"].write(rp.dispHeaderFields(all_fields)+"\t"+header_recompute+"\n")

    names = None
    if data is not None and data.hasNames() and "queries_named" in filenames:
        names = data.getNames()
        filesfp["queries_named"] = open(filenames["queries_named"], mode)
        if with_headers:
            filesfp["queries_named"].write(rp.dispHeaderFields(all_fields)+"\t"+header_recompute+"\n")
    
    if "support" in filenames:
        filesfp["support"] = open(filenames["support"], mode)
        
    #### TO DEBUG: output all shown in siren, i.e. no filtering
    ## for pos in range(len(results_batch["batch"])):
    addto = ""
    for pos in results_batch["results"]:
        org = results_batch["batch"][pos]
        
        if data_recompute is not None:
            red = org.copy()
            red.recompute(data_recompute)
            acc_diff = (red.getAcc()-org.getAcc())/org.getAcc()
            addto = "\t"+red.disp(list_fields=fields_recompute)+"\t%f" % acc_diff
        filesfp["queries"].write(org.disp(list_fields=all_fields)+addto+'\n')
        if filesfp["queries_named"] is not None:
            filesfp["queries_named"].write(org.disp(names, list_fields=all_fields_named)+addto+'\n')
        if filesfp["support"] is not None:
            filesfp["support"].write(org.dispSupp()+'\n')

    for (ffi, ffp) in filesfp.items():
        if ffp is not None and filenames.get(ffi, "") != "-":
            ffp.close()

def loadPackage(filename, pm):

    package = Package(filename)
    elements_read = package.read(pm)        

    if elements_read.get("data") is not None:
        data = elements_read.get("data")
    else:
        data = None
    if elements_read.get("preferences"):
        params = elements_read.get("preferences")
    else:
        params = None

    return params, data

def applyVarsMask(data, params):
    params_l = turnToDict(params)
    return data.applyDisableMasks(params_l.get("mask_vars_LHS"), params_l.get("mask_vars_RHS"), params_l.get("mask_rows"))

def run(args):
    
    loaded = loadAll(args)
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"]) 
    miner = instMiner(data, params, logger)
    try:
        miner.full_run()
    except KeyboardInterrupt:
        miner.initial_pairs.saveToFile()
        logger.printL(1, 'Stopped...', "log")

    outputResults(filenames, miner.final, data)
    logger.clockTac(0, None)

###############################
 # def run_filterRM(args):
    
 #    loaded = loadAll(args)
 #    params, data, logger, filenames, reds = (loaded["params"], loaded["data"], loaded["logger"],
 #                                             loaded["filenames"], loaded["reds"]) 

 #    constraints = Constraints(data, params)

 #    candidate_ids = range(len(reds))
 #    scores = numpy.zeros((len(reds)+2, len(reds)))
 #    keep_ids = []
 #    rm = RedModel(data)
 #    best = (0, -1)
 #    while best[-1] is not None:
 #        best = (0, None)
 #        tic = datetime.datetime.now()
 #        for rr, ri in enumerate(candidate_ids):
 #            top = rm.getTopDeltaRed(reds[ri], data)
 #            scores[ri, len(keep_ids)] = top[0]
 #            if top[0] < best[0]:
 #                best = (top[0], top[1], rr)
 #                # print top, reds[ri].disp()

 #        if best[-1] is not None:
 #            ri = candidate_ids.pop(best[-1])
 #            scores[-2, len(keep_ids)] = (datetime.datetime.now()-tic).total_seconds()
 #            scores[-1, len(keep_ids)] = ri
 #            keep_ids.append(ri)
 #            rm.addRed(reds[ri], data, best[1])
 #            print "%f\t%d\t%s" % (best[0], ri, reds[ri].disp())

 #    numpy.savetxt('scores.txt', scores, fmt="%f")
###############################

def run_filter(args):

    loaded = loadAll(args)
    params, data, logger, filenames, reds = (loaded["params"], loaded["data"], loaded["logger"],
                                             loaded["filenames"], loaded["reds"]) 

    constraints = Constraints(data, params)
    rp = Redescription.getRP()
    reds = []
    with open("/home/galbrun/current/redescriptions.csv") as fd:
        rp.parseRedList(fd, data, reds)

    rr_tests = [[1, 32, 6, 5, 29, 94], [23, 12], [7, 66, 11, 29]]
    rr_tests = [[73]] ## [2]
    for ri, add_redids in enumerate(rr_tests):

        include_redids = [84, 77, 53, 29, 94]

        bbatch = Batch([reds[i] for i in include_redids]+[reds[i] for i in add_redids])
        org_ids = bbatch.selected(constraints.getActions("final"))
        
        batch = Batch([reds[i] for i in include_redids])
        pids = batch.selected(constraints.getActions("final"))
        batch.extend([reds[i] for i in add_redids])
        # tmp_ids = batch.selected(self.constraints.getActions("redundant"))
        ticc = datetime.datetime.now()
        new_ids = range(len(include_redids), len(include_redids)+len(add_redids))
        tmp_ids = batch.selected(constraints.getActions("final"), ids= pids+new_ids, new_ids=new_ids)
        tacc = datetime.datetime.now()
        print "Elapsed ", ri, tacc-ticc
        if tmp_ids != org_ids:
            print "Not identical"
        pdb.set_trace()
        print len(tmp_ids), len(org_ids)
        
    return [batch[i] for i in tmp_ids]

    ## miner = instMiner(data, params, logger)
    ## try:
    ##     miner.full_run()
    ## except KeyboardInterrupt:
    ##     logger.printL(1, 'Stopped...', "log")
    ## 
    ## outputResults(filenames, miner.final, data)
    ## logger.clockTac(0, None)


def run_splits(args, splt=""):
    nb_splits = 5
    tmp = re.match("splits(?P<nbs>[0-9]+)\s*", splt)
    if tmp is not None:
        nb_splits = int(tmp.group("nbs"))
        
    loaded = loadAll(args)
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"]) 
    if "package" in filenames:
        parts = filenames["package"].split("/")[-1].split(".")
        pp = filenames["basis"].split("/")
        pp[-1] = ".".join(parts[:-1])
        filenames["basis"] = "/".join(pp)
    fold_cols = data.findCandsFolds(strict=True)

    if len(fold_cols) == 0:
        fold_cols = [None]
    else:
        for fci in fold_cols:
            data.cols[fci[0]][fci[1]].setDisabled()

    for fci in fold_cols:
        if fci is None:
            logger.printL(2, "Data has no folds, generating...", "log")
            sss = data.getSplit(nbsubs=nb_splits)
            data.addFoldsCol()
            suff = "rand"
            splt_pckgf = filenames["basis"]+ ("_split-%d:%s_empty.siren" % (nb_splits, suff))
            saveAsPackage(splt_pckgf, data, preferences=params, pm=loaded["pm"])        
        else:
            logger.printL(2, "Using existing fold: side %s col %s" % fci, "log")
            sss = data.extractFolds(fci[0], fci[1])
            nb_splits = len(sss)
            suff = data.cols[fci[0]][fci[1]].getName()
        print "SIDS", suff, sorted(data.getFoldsInfo()["split_ids"].items(), key=lambda x: x[1])
        print data
        splt_pckgf = filenames["basis"]+ ("_split-%d:%s.siren" % (nb_splits, suff))
        splt_statf = filenames["basis"]+ ("_split-%d:%s.txt" % (nb_splits, suff))            

        stM = StatsMiner(data, params, logger)
        reds_list, all_stats, summaries, list_fields, stats_fields = stM.run_stats()
        
        splt_fk = filenames["basis"]+ ("_split-%d:%s-kall.txt" % (nb_splits, suff))            
        with open(splt_fk, "w") as f:
            f.write(rp.printRedList(reds_list, fields=list_fields+["track"]))
            
        for fk, dt in summaries.items():
            splt_fk = filenames["basis"]+ ("_split-%d:%s-k%d.txt" % (nb_splits, suff, fk))            
            with open(splt_fk, "w") as f:
                f.write(rp.printRedList(dt["reds"], fields=list_fields+["track"]))
                
        nbreds = numpy.array([len(ll) for (li, ll) in all_stats.items() if li > -1])
        tot = numpy.array(all_stats[-1])
        if nbreds.sum() > 0:
            summary_mat = numpy.hstack([numpy.vstack([tot.min(axis=0), tot.max(axis=0), tot.mean(axis=0), tot.std(axis=0)]), numpy.array([[nbreds.min()], [nbreds.max()], [nbreds.mean()], [nbreds.std()]])])

            info_plus = "\nrows:min\tmax\tmean\tstd\tnb_folds:%d" % (len(all_stats)-1)
            numpy.savetxt(splt_statf, summary_mat, fmt="%f", delimiter="\t", header="\t".join(stats_fields+["nb reds"])+info_plus)
            # saveAsPackage(splt_pckgf, data, preferences=params, pm=loaded["pm"], reds=reds_list)
        else:
            with open(splt_statf, "w") as fo:
                fo.write("No redescriptions found")
        # for red in reds_list:
        #     print red.disp()


def run_printout(args):
    suff = args[-1].strip("printout")
    if len(suff) == 0:
        suff = "_reprint"
    loaded = loadAll(args[:-1])
    params, data, logger, filenames, reds = (loaded["params"], loaded["data"], loaded["logger"],
                                             loaded["filenames"], loaded["reds"])
    rp = Redescription.getRP()
    qfilename = None
    if reds is None and "queries" in filenames:
        qfilename = filenames["queries"]        
    if "queries_second" in filenames:
        qfilename = filenames["queries_second"]

    if qfilename is not None:
        reds = []
        try:
            with open(qfilename) as fd:
                rp.parseRedList(fd, data, reds)
        except IOError:
            reds = []
    
    #### OUT
    parts = filenames["queries"].split(".")
    if len(parts) > 1:
        if "." in suff:
            filename = ".".join(parts[:-2] + [parts[-2]+ suff])
        else:
            filename = ".".join(parts[:-2] + [parts[-2]+ suff, parts[-1]])
    else:
        filename = filenames["queries"] + suff
    print "FILENAME", filename
    if type(reds) is list and len(reds) > 0 and type(reds[0]) is dict and "items" in reds[0]:
        red_contents = []
        for r in reds:
            red_contents.extend([r["items"][rid] for rid in r["rshowids"]])
    else:
        red_contents = reds
                
    params = getPrintParams(filename, data)
    writeRedescriptions(red_contents, filename, **params)
                

def run_rnd(args):

    pref_dir = os.path.dirname(os.path.abspath(__file__))
    conf_defs = [pref_dir + "/miner_confdef.xml", pref_dir + "/inout_confdef.xml", pref_dir + "/rnd_confdef.xml"]
    
    loaded = loadAll(args, conf_defs)
    params, data, logger, filenames = (loaded["params"], loaded["data"], loaded["logger"], loaded["filenames"])

    params_l = turnToDict(params)
    select_red = None 
    if len(params_l.get("select_red", "")) > 0:
        select_red = params_l["select_red"]
    prec_all = None
    if params_l.get("agg_prec", -1) >= 0:
        prec_all = params_l["agg_prec"]
    count_vname = params_l.get("count_vname", "COUNTS")
            
    rf = RndFactory(org_data=data)    
    with_traits=False
    if "traits_data" in filenames:
        traits_data = Data([filenames["traits_data"], None]+filenames["add_info"], filenames["style_data"])
        rf.setTraits(traits_data)
        with_traits=True
        
    if params_l.get("rnd_seed", -1) >= 0:
        rf.setSeed(params_l["rnd_seed"])

    stop = False
    for rnd_meth in params_l["rnd_meth"]:
        nb_copies = params_l["rnd_series_size"]
        if rnd_meth == "none":
            nb_copies = 1
        
        for i in range(nb_copies):
            sub_filenames = dict(filenames)
            suff = "_%s-%d" % (rnd_meth, i)
            sub_filenames["basis"] += suff
            for k in ["queries", "queries_named", "support"]:
                if k in sub_filenames:
                    parts = sub_filenames[k].split(".")
                    parts[-2] += suff
                    sub_filenames[k] = ".".join(parts)

            Dsub, sids, back, store = rf.makeupRndData(rnd_meth=rnd_meth, with_traits=with_traits, count_vname=count_vname, select_red=select_red, prec_all=prec_all)
            logger.printL(2, "STARTING Random series %s %d" % (rnd_meth, i), "log")
            logger.printL(2, Dsub, "log")

            miner = instMiner(Dsub, params, logger)
            try:
                miner.full_run()
            except KeyboardInterrupt:
                miner.initial_pairs.saveToFile()
                logger.printL(1, 'Stopped...', "log")
                stop = True
                
            outputResults(sub_filenames, miner.final, Dsub)
            logger.clockTac(0, None)
            if stop:
                exit()
    
##### MAIN
###########
    
if __name__ == "__main__":

    if re.match("printout", sys.argv[-1]):
        run_printout(sys.argv)
    elif re.match("rnd", sys.argv[-1]):
        run_rnd(sys.argv[:-1])
    elif re.match("splits", sys.argv[-1]):
        run_splits(sys.argv[:-1], sys.argv[-1])
    elif sys.argv[-1] == "filter":
        run_filter(sys.argv[:-1])
    elif sys.argv[-1] == "filterRM":
        run_filterRM(sys.argv[:-1])
    else:
        run(sys.argv)
