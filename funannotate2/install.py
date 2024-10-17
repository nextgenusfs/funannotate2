import os
import io
import sys
import hashlib
import json
import subprocess
import datetime
import requests
import pyhmmer
from urllib.request import urlopen
import xml.etree.cElementTree as cElementTree
from .log import startLogging, system_info, finishLogging
from .config import env
from .utilities import download, load_json, runSubprocess, checkfile
from .fastx import countfasta


def install(args):
    # start logger
    global logger
    logger = startLogging()
    log = logger.info
    system_info(log)

    # now get databases, etc
    if not os.path.isdir(env["FUNANNOTATE2_DB"]):
        os.makedirs(env["FUNANNOTATE2_DB"])
    logger.info(
        f"The backend database location is from the $FUNANNOTATE2_DB env variable: {env['FUNANNOTATE2_DB']}"
    )

    global today
    today = datetime.datetime.today().strftime("%Y-%m-%d")

    db2install = []
    if "all" in args.db:
        db2install = [
            "merops",
            "uniprot",
            "dbCAN",
            "pfam",
            "go",
            "mibig",
            "interpro",
            "gene2product",
        ]
    else:
        db2install = args.db

    # load download links from gitlab if possible, need to change this when merge into master
    global DBURL
    try:
        logger.info("Retrieving download links from GitHub Repo")
        DBURL = json.loads(
            requests.get(
                "https://raw.githubusercontent.com/nextgenusfs/funannotate2/master/funannotate2/downloads.json"
            ).text
        )["downloads"]
    except:  # noqa: E722
        logger.error(
            "Unable to download links from GitHub, using local copy from funannotate2"
        )
        DBURL = load_json(os.path.join(os.path.dirname(__file__), "downloads.json"))[
            "downloads"
        ]
    logger.info(json.dumps(DBURL, indent=2))

    # pull the database JSON file
    # if text file with DB info is in database folder, parse into Dictionary
    DatabaseFile = os.path.join(env["FUNANNOTATE2_DB"], "funannotate-db-info.json")
    dbinfo = {}
    if checkfile(DatabaseFile):
        dbinfo = load_json(DatabaseFile)

    # let user know what you are doing
    if args.force is False:
        if args.update is True:
            logger.info("Checking for newer versions of database files")
    else:
        logger.info(
            "--force option is passed, will re-download all database and re-format"
        )

    # now go through installdbs and download/update/etc
    for x in db2install:
        if x == "uniprot":
            data = {}
            if "uniprot" in dbinfo:
                data = dbinfo.get("uniprot")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("uniprot"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = uniprotDB(wget=args.wget)
                dbinfo["uniprot"] = data

        elif x == "merops":
            data = {}
            if "merops" in dbinfo:
                data = dbinfo.get("merops")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("merops"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = meropsDB(wget=args.wget)
                dbinfo["merops"] = data

        elif x == "dbCAN":
            data = {}
            if "dbCAN" in dbinfo:
                data = dbinfo.get("dbCAN")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("dbCAN"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = dbCANDB(wget=args.wget)
                dbinfo["dbCAN"] = data

        elif x == "pfam":
            data = {}
            if "pfam" in dbinfo:
                data = dbinfo.get("pfam")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("pfam"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = pfamDB(wget=args.wget)
                dbinfo["pfam"] = data

        elif x == "go":
            data = {}
            if "go" in dbinfo:
                data = dbinfo.get("go")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("go"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = goDB(wget=args.wget)
                dbinfo["go"] = data

        elif x == "interpro":
            data = {}
            if "interpro" in dbinfo:
                data = dbinfo.get("interpro")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("interpro"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = interproDB(wget=args.wget)
                dbinfo["interpro"] = data

        elif x == "mibig":
            data = {}
            if "mibig" in dbinfo:
                data = dbinfo.get("mibig")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("mibig"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = mibigDB(wget=args.wget)
                dbinfo["mibig"] = data

        elif x == "gene2product":
            data = {}
            if "gene2product" in dbinfo:
                data = dbinfo.get("gene2product")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                if remote_md5_updated(DBURL.get("gene2product"), data["md5"]):
                    data = {}
            # if its empty then need to run it
            if not data:
                data = curatedDB(wget=args.wget)
                dbinfo["gene2product"] = data

    # for now just output database file as json
    logger.info(f"Database installation complete:\n{json.dumps(dbinfo, indent=2)}")

    # write the database json file
    with open(DatabaseFile, "w") as dbfile:
        json.dump(dbinfo, dbfile, indent=2)

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])


def calcmd5(file):
    md5local = None
    with open(file, "rb") as infile:
        data = infile.read()
        md5local = hashlib.md5(data).hexdigest()
    return md5local


def calcmd5remote(url, max_file_size=100 * 1024 * 1024):
    remote = urlopen(url)
    hash = hashlib.md5()
    total_read = 0
    while True:
        data = remote.read(4096)
        total_read += 4096
        if not data or total_read > max_file_size:
            break
        hash.update(data)
    return hash.hexdigest()


def remote_md5_updated(url, md5):
    newmd5 = calcmd5remote(url)
    if newmd5 == md5:
        return False
    else:
        return True


def meropsDB(wget=False):
    fasta = os.path.join(env["FUNANNOTATE2_DB"], "meropsscan.lib")
    filtered = os.path.join(env["FUNANNOTATE2_DB"], "merops.formatted.fa")
    database = os.path.join(env["FUNANNOTATE2_DB"], "merops.dmnd")
    logger.info("Downloading Merops database")
    for x in [fasta, filtered, database]:
        if checkfile(x):
            os.remove(x)
    # download here
    download(DBURL.get("merops"), fasta, wget=wget)
    # calc md5 checksum of raw input/download
    md5 = calcmd5(fasta)
    # reformat fasta headers
    with open(filtered, "w") as filtout:
        with io.open(fasta, encoding="utf8", errors="ignore") as infile:
            for line in infile:
                if line.startswith(">"):
                    line = line.rstrip()
                    ID = line.split()[0]
                    family = line.split("#")[1]
                    filtout.write("{:} {:}\n".format(ID, family))
                else:
                    filtout.write(line)
    cmd = ["diamond", "makedb", "--in", "merops.formatted.fa", "--db", "merops"]
    logger.info(f"Building diamond database: {' '.join(cmd)}")
    runSubprocess(cmd, os.path.join(env["FUNANNOTATE2_DB"]), logger)
    num_records = countfasta(filtered)
    return {
        "type": "diamond",
        "db": database,
        "version": "12.5",
        "date": "2023-01-19",
        "n_records": num_records,
        "md5": md5,
    }


def uniprotDB(wget=False):
    fasta = os.path.join(env["FUNANNOTATE2_DB"], "uniprot_sprot.fasta")
    database = os.path.join(env["FUNANNOTATE2_DB"], "uniprot.dmnd")
    versionfile = os.path.join(env["FUNANNOTATE2_DB"], "uniprot.release-date.txt")
    # remove existing files if they exist
    logger.info("Downloading UniProtKB/Swiss-Prot database")
    for x in [fasta, fasta + ".gz", versionfile, database]:
        if checkfile(x):
            os.remove(x)
    download(DBURL.get("uniprot"), fasta + ".gz", wget=wget)
    download(DBURL.get("uniprot-release"), versionfile, wget=wget)
    subprocess.call(
        ["gunzip", "-f", "uniprot_sprot.fasta.gz"],
        cwd=os.path.join(env["FUNANNOTATE2_DB"]),
    )
    md5 = calcmd5(versionfile)
    unidate = ""
    univers = ""
    with io.open(versionfile, encoding="utf8", errors="ignore") as infile:
        for line in infile:
            if line.startswith("UniProtKB/Swiss-Prot Release"):
                rest, datepart = line.split(" of ")
                unidate = datetime.datetime.strptime(
                    datepart.rstrip(), "%d-%b-%Y"
                ).strftime("%Y-%m-%d")
                univers = rest.split(" ")[-1]
    # build diamond database
    cmd = ["diamond", "makedb", "--in", "uniprot_sprot.fasta", "--db", "uniprot"]
    logger.info(f"Building diamond database: {' '.join(cmd)}")
    runSubprocess(cmd, os.path.join(env["FUNANNOTATE2_DB"]), logger)
    num_records = countfasta(fasta)
    return {
        "type": "diamond",
        "db": database,
        "version": univers,
        "date": unidate,
        "n_records": num_records,
        "md5": md5,
    }


def dbCANDB(wget=False):
    hmm = os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.hmm")
    familyinfo = os.path.join(env["FUNANNOTATE2_DB"], "dbCAN-fam-HMMs.txt")
    versionfile = os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.changelog.txt")
    logger.info("Downloading dbCAN database")
    # delete existint
    for x in [
        os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp"),
        hmm,
        familyinfo,
        versionfile,
    ]:
        if checkfile(x):
            os.remove(x)
    # download data
    download(
        DBURL.get("dbCAN"), os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp"), wget=wget
    )
    download(DBURL.get("dbCAN-tsv"), familyinfo, wget=wget)
    download(DBURL.get("dbCAN-log"), versionfile, wget=wget)
    md5 = calcmd5(os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp"))
    num_records = 0
    dbdate = ""
    dbvers = ""
    with open(hmm, "w") as out:
        with io.open(
            os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp"),
            encoding="utf8",
            errors="ignore",
        ) as input:
            for line in input:
                if line.startswith("NAME"):
                    num_records += 1
                    line = line.replace(".hmm\n", "\n")
                out.write(line)
    with io.open(versionfile, encoding="utf8", errors="ignore") as infile:
        head = [next(infile) for x in range(2)]
    dbdate = head[1].replace("# ", "").rstrip()
    dbvers = head[0].split(" ")[-1].rstrip()
    dbdate = datetime.datetime.strptime(dbdate, "%m/%d/%Y").strftime("%Y-%m-%d")
    logger.info("Creating dbCAN HMM database and pressing with pyhmmer")
    pyhmmer.hmmer.hmmpress(pyhmmer.plan7.HMMFile(hmm), hmm)
    if checkfile(os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp")):
        os.remove(os.path.join(env["FUNANNOTATE2_DB"], "dbCAN.tmp"))
    return {
        "type": "hmmer3",
        "db": hmm,
        "version": dbvers,
        "date": dbdate,
        "n_records": num_records,
        "md5": md5,
    }


def pfamDB(wget=False):
    hmm = os.path.join(env["FUNANNOTATE2_DB"], "Pfam-A.hmm")
    familyinfo = os.path.join(env["FUNANNOTATE2_DB"], "Pfam-A.clans.tsv")
    versionfile = os.path.join(env["FUNANNOTATE2_DB"], "Pfam.version")
    # delete previous if exists
    for x in [
        hmm,
        hmm + ".gz",
        familyinfo,
        familyinfo + ".gz",
        versionfile,
        versionfile + ".gz",
    ]:
        if checkfile(x):
            os.remove(x)
    logger.info("Downloading Pfam database")
    download(DBURL.get("pfam"), hmm + ".gz", wget=wget)
    download(DBURL.get("pfam-tsv"), familyinfo + ".gz", wget=wget)
    download(DBURL.get("pfam-log"), versionfile + ".gz", wget=wget)
    subprocess.call(
        ["gunzip", "-f", "Pfam-A.hmm.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"])
    )
    subprocess.call(
        ["gunzip", "-f", "Pfam-A.clans.tsv.gz"],
        cwd=os.path.join(env["FUNANNOTATE2_DB"]),
    )
    md5 = calcmd5(versionfile + ".gz")
    subprocess.call(
        ["gunzip", "-f", "Pfam.version.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"])
    )
    num_records = 0
    pfamdate = ""
    pfamvers = ""
    with io.open(versionfile, encoding="utf8", errors="ignore") as input:
        for line in input:
            if line.startswith("Pfam release"):
                pfamvers = line.split(": ")[-1].rstrip()
            if line.startswith("Pfam-A families"):
                num_records = int(line.split(": ")[-1].rstrip())
            if line.startswith("Date"):
                pfamdate = line.split(": ")[-1].rstrip()
    logger.info("Creating Pfam HMM database and pressing with pyhmmer")
    pyhmmer.hmmer.hmmpress(pyhmmer.plan7.HMMFile(hmm), hmm)
    return {
        "type": "hmmer3",
        "db": hmm,
        "version": pfamvers,
        "date": pfamdate,
        "n_records": num_records,
        "md5": md5,
    }


def goDB(wget=False):
    goOBO = os.path.join(env["FUNANNOTATE2_DB"], "go.obo")
    logger.info("Downloading GO Ontology database")
    for x in [goOBO]:
        if checkfile(x):
            os.remove(x)
    download(DBURL.get("go"), goOBO, wget=wget)
    md5 = calcmd5(goOBO)
    num_records = 0
    version = ""
    with io.open(goOBO, encoding="utf8", errors="ignore") as infile:
        for line in infile:
            if line.startswith("data-version:"):
                version = line.split(" ")[1].rstrip().replace("releases/", "")
            if line.startswith("[Term]"):
                num_records += 1
    return {
        "type": "text",
        "db": goOBO,
        "version": version,
        "date": version,
        "n_records": num_records,
        "md5": md5,
    }


def mibigDB(wget=False):
    fasta = os.path.join(env["FUNANNOTATE2_DB"], "mibig.fa")
    database = os.path.join(env["FUNANNOTATE2_DB"], "mibig.dmnd")
    logger.info("Downloading MiBIG Secondary Metabolism database")
    for x in [fasta, database]:
        if checkfile(x):
            os.remove(x)
    download(DBURL.get("mibig"), fasta, wget=wget)
    md5 = calcmd5(fasta)
    version = os.path.basename(DBURL.get("mibig")).split("_")[-1].replace(".fasta", "")
    cmd = ["diamond", "makedb", "--in", "mibig.fa", "--db", "mibig"]
    logger.info(f"Building diamond database: {' '.join(cmd)}")
    runSubprocess(cmd, os.path.join(env["FUNANNOTATE2_DB"]), logger)
    num_records = countfasta(fasta)
    return {
        "type": "diamond",
        "db": database,
        "version": version,
        "date": today,
        "n_records": num_records,
        "md5": md5,
    }


def interproDB(wget=False):
    iprXML = os.path.join(env["FUNANNOTATE2_DB"], "interpro.xml")
    iprTSV = os.path.join(env["FUNANNOTATE2_DB"], "interpro.tsv")
    logger.info("Downloading InterProScan Mapping file")
    for x in [iprXML, iprTSV, iprXML + ".gz"]:
        if checkfile(x):
            os.remove(x)
    download(DBURL.get("interpro"), iprXML + ".gz", wget=wget)
    download(DBURL.get("interpro-tsv"), iprTSV, wget=wget)
    md5 = calcmd5(iprXML + ".gz")
    subprocess.call(
        ["gunzip", "-f", "interpro.xml.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"])
    )
    num_records = ""
    version = ""
    iprdate = ""
    for event, elem in cElementTree.iterparse(iprXML):
        if elem.tag == "release":
            for x in list(elem):
                if x.attrib["dbname"] == "INTERPRO":
                    num_records = int(x.attrib["entry_count"])
                    version = x.attrib["version"]
                    iprdate = x.attrib["file_date"]
    try:
        iprdate = datetime.datetime.strptime(iprdate, "%d-%b-%y").strftime("%Y-%m-%d")
    except ValueError:
        iprdate = datetime.datetime.strptime(iprdate, "%d-%b-%Y").strftime("%Y-%m-%d")
    return {
        "type": "xml",
        "db": iprXML,
        "version": version,
        "date": iprdate,
        "n_records": num_records,
        "md5": md5,
    }


def curatedDB(wget=False):
    curatedFile = os.path.join(env["FUNANNOTATE2_DB"], "ncbi_cleaned_gene_products.txt")
    logger.info("Downloaded curated gene names and product descriptions")
    for x in [curatedFile]:
        if checkfile(x):
            os.remove(x)
    download(DBURL.get("gene2product"), curatedFile, wget=wget)
    md5 = calcmd5(curatedFile)
    num_records = 0
    curdate = ""
    version = ""
    with io.open(curatedFile, encoding="utf8", errors="ignore") as infile:
        for line in infile:
            if line.startswith("#version"):
                version = line.split(" ")[-1].rstrip()
            elif line.startswith("#Date"):
                curdate = line.split(" ")[-1].rstrip()
            else:
                num_records += 1
    curdate = datetime.datetime.strptime(curdate, "%m-%d-%Y").strftime("%Y-%m-%d")
    return {
        "type": "text",
        "db": curatedFile,
        "version": version,
        "date": today,
        "n_records": num_records,
        "md5": md5,
    }
