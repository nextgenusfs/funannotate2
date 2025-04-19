import datetime
import hashlib
import io
import json
import os
import subprocess
import sys
import xml.etree.cElementTree as cElementTree
from urllib.request import urlopen

import pyhmmer
import requests

from .config import env
from .fastx import countfasta
from .log import finishLogging, startLogging, system_info
from .utilities import checkfile, download, load_json


def install(args):
    """
    Install databases based on user input and update options.

    This function initializes logging and retrieves database download links. It checks
    for existing databases and updates them if necessary, based on user-specified options
    such as force installation or update checks. The function supports installing various
    databases like UniProt, MEROPS, dbCAN, Pfam, and others, and logs the installation
    process and results.

    If the --show option is provided, it displays the currently installed databases
    without performing any installations or updates.

    Parameters:
    - args: Command-line arguments provided by the user.

    Returns:
    - None
    """
    # start logger
    global logger
    logger = startLogging()
    log = logger.info
    system_info(log)

    # now get databases and ensure folders exist
    if not os.path.isdir(env["FUNANNOTATE2_DB"]):
        os.makedirs(env["FUNANNOTATE2_DB"])
    logger.info(
        f"The backend database location is from the $FUNANNOTATE2_DB env variable: {env['FUNANNOTATE2_DB']}"
    )
    if not os.path.isdir(os.path.join(env["FUNANNOTATE2_DB"], "pretrained")):
        os.makedirs(os.path.join(env["FUNANNOTATE2_DB"], "pretrained"))

    global today
    today = datetime.datetime.today().strftime("%Y-%m-%d")

    # pull the database JSON file
    # if text file with DB info is in database folder, parse into Dictionary
    DatabaseFile = os.path.join(env["FUNANNOTATE2_DB"], "funannotate-db-info.json")
    dbinfo = {}
    if checkfile(DatabaseFile):
        dbinfo = load_json(DatabaseFile)

    # If --show option is provided, display the installed databases as JSON and exit
    if hasattr(args, "show") and args.show:
        if not dbinfo:
            logger.info("No databases are currently installed.")
        else:
            # Clean up the database paths to show only filenames for cleaner output
            display_info = {}
            for db_name, db_info in dbinfo.items():
                display_info[db_name] = db_info.copy()
                if "db" in display_info[db_name]:
                    display_info[db_name]["db"] = os.path.basename(display_info[db_name]["db"])

            # Print as formatted JSON
            logger.info("Currently installed databases:")
            logger.info(json.dumps(display_info, indent=2, sort_keys=True))
        # Finish logging and exit
        finishLogging(log, vars(sys.modules[__name__])["__name__"])
        return

    # Check if --db argument is provided for installation
    if not hasattr(args, "db") or not args.db:
        logger.error(
            "No databases specified for installation. Use -d/--db to specify databases to install."
        )
        finishLogging(log, vars(sys.modules[__name__])["__name__"])
        return

    # Proceed with installation
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
            "mito",
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
        logger.error("Unable to download links from GitHub, using local copy from funannotate2")
        DBURL = load_json(os.path.join(os.path.dirname(__file__), "downloads.json"))["downloads"]
    logger.info(json.dumps(DBURL, indent=2))

    # let user know what you are doing
    if args.force is False:
        if args.update is True:
            logger.info("Checking for newer versions of database files")
    else:
        logger.info("--force option is passed, will re-download all database and re-format")

    # now go through installdbs and download/update/etc
    for x in db2install:
        if x == "uniprot":
            data = {}
            if "uniprot" in dbinfo:
                data = dbinfo.get("uniprot")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                check = remote_md5_updated(DBURL.get("uniprot-release"), data["md5"])
                if check is not False:
                    logger.info(
                        f"UniProtKB checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in UniProtKB/Swiss-Prot database")
            # if its empty then need to run it
            if not data:
                data = uniprotDB(wget=args.wget)
                dbinfo["uniprot"] = data

        elif x == "mito":
            data = {}
            if "mito" in dbinfo:
                data = dbinfo.get("mito")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                check = remote_md5_updated(DBURL.get("mito-release"), data["md5"])
                if check is not False:
                    logger.info(
                        f"UniProtKB checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in UniProtKB/Swiss-Prot database")
            # if its empty then need to run it
            if not data:
                data = mitoDB(wget=args.wget)
                dbinfo["mito"] = data

        elif x == "merops":
            data = {}
            if "merops" in dbinfo:
                data = dbinfo.get("merops")
            if args.force is True:
                data = {}
            elif args.update is True and "md5" in data:
                check = remote_md5_updated(DBURL.get("merops"), data["md5"])
                if check is not False:
                    logger.info(
                        f"MEROPS checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in MEROPS protease database")
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
                check = remote_md5_updated(DBURL.get("dbCAN"), data["md5"])
                if check is not False:
                    logger.info(
                        f"dbCAN checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in dbCAN database")
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
                check = remote_md5_updated(DBURL.get("pfam-log"), data["md5"])
                if check is not False:
                    logger.info(
                        f"Pfam-A checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in Pfam-A database")
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
                check = remote_md5_updated(DBURL.get("go"), data["md5"])
                if check is not False:
                    logger.info(
                        f"GO-OBO checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in Gene Ontology database")
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
                check = remote_md5_updated(DBURL.get("interpro-tsv"), data["md5"])
                if check is not False:
                    logger.info(
                        f"InterPro checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in InterPro database")
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
                check = remote_md5_updated(DBURL.get("mibig"), data["md5"])
                if check is not False:
                    logger.info(
                        f"MiBig checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in MiBiG database")
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
                check = remote_md5_updated(DBURL.get("gene2product"), data["md5"])
                if check is not False:
                    logger.info(
                        f"Gene2Product checksum suggests an update: existing={data['md5']} remote={check}"
                    )
                    data = {}
                else:
                    logger.info("No change detected in Gene2Product database")
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
    """
    Calculate the MD5 hash value of a file.

    This function reads the contents of a specified file in binary mode and calculates
    its MD5 hash value. The hash is returned as a hexadecimal string.

    Parameters:
    - file (str): The path to the file to calculate the MD5 hash for.

    Returns:
    - str: The MD5 hash value of the file.
    """
    md5local = None
    with open(file, "rb") as infile:
        data = infile.read()
        md5local = hashlib.md5(data).hexdigest()
    return md5local


def calcmd5remote_original(url, max_file_size=100 * 1024 * 1024):
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


def calcmd5remote(url, timeout=60):
    """
    Calculate the MD5 hash of a remote file's headers or a small sample.

    This function tries HTTPS first, then falls back to FTP if needed.

    Parameters:
    - url (str): The URL of the file to calculate the MD5 hash for.
    - timeout (int, optional): Timeout in seconds for the request (default is 60).

    Returns:
    - str: The MD5 hash value or a header-based identifier of the file.
    """
    # Try HTTPS first
    try:
        # First try to get headers only, which is much faster
        response = requests.head(url, timeout=timeout, verify=False)
        response.raise_for_status()

        # Check if we have ETag or Last-Modified headers
        etag = response.headers.get("ETag")
        if etag:
            # Remove quotes if present
            return etag.strip("\"'")

        last_modified = response.headers.get("Last-Modified")
        if last_modified:
            # Use Last-Modified as a proxy for content
            return hashlib.md5(last_modified.encode()).hexdigest()

        # If we don't have useful headers, fall back to sampling the file
        response = requests.get(url, stream=True, timeout=timeout, verify=False)
        response.raise_for_status()

        hash = hashlib.md5()
        total_read = 0
        max_sample = 1024 * 1024  # 1MB sample

        for chunk in response.iter_content(chunk_size=4096):
            if not chunk:
                break
            total_read += len(chunk)
            hash.update(chunk)
            if total_read >= max_sample:
                break

        return hash.hexdigest()

    except requests.exceptions.RequestException as e:
        print(f"HTTPS request failed: {str(e)}")

        # If HTTPS fails and the URL is using HTTPS to access an FTP server, try direct FTP
        if "ftp." in url and url.startswith("https://"):
            try:
                # Convert HTTPS URL to FTP URL
                ftp_url = url.replace("https://", "ftp://")
                print(f"Trying FTP fallback: {ftp_url}")

                import socket
                from urllib.request import urlopen

                # Set a socket timeout for FTP connections
                socket.setdefaulttimeout(timeout)

                remote = urlopen(ftp_url)
                hash = hashlib.md5()
                total_read = 0
                max_sample = 1024 * 1024  # 1MB sample

                while True:
                    data = remote.read(4096)
                    if not data or total_read >= max_sample:
                        break
                    total_read += len(data)
                    hash.update(data)

                return hash.hexdigest()
            except Exception as ftp_e:
                print(f"FTP fallback also failed: {str(ftp_e)}")

        # Return a placeholder value that won't match any real MD5
        return "ERROR_CALCULATING_MD5"


def remote_md5_updated(url, md5):
    newmd5 = calcmd5remote(url)
    if newmd5 == md5:
        return False
    else:
        return newmd5


def meropsDB(wget=False):
    fasta = os.path.join(env["FUNANNOTATE2_DB"], "meropsscan.lib")
    filtered = os.path.join(env["FUNANNOTATE2_DB"], "merops.formatted.fa")
    database = os.path.join(env["FUNANNOTATE2_DB"], "merops.dmnd")
    logger.info("Downloading Merops database")
    for x in [fasta, filtered, database]:
        if checkfile(x):
            os.remove(x)
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("merops"))
    # download here
    download(DBURL.get("merops"), fasta, wget=wget)
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
    subprocess.call(cmd, cwd=os.path.join(env["FUNANNOTATE2_DB"]))
    # runSubprocess(cmd, os.path.join(env['FUNANNOTATE2_DB']), logger)
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
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("uniprot-release"))
    download(DBURL.get("uniprot"), fasta + ".gz", wget=wget)
    download(DBURL.get("uniprot-release"), versionfile, wget=wget)
    subprocess.call(
        ["gunzip", "-f", "uniprot_sprot.fasta.gz"],
        cwd=os.path.join(env["FUNANNOTATE2_DB"]),
    )
    unidate = ""
    univers = ""
    with io.open(versionfile, encoding="utf8", errors="ignore") as infile:
        for line in infile:
            if line.startswith("UniProtKB/Swiss-Prot Release"):
                rest, datepart = line.split(" of ")
                unidate = datetime.datetime.strptime(datepart.rstrip(), "%d-%b-%Y").strftime(
                    "%Y-%m-%d"
                )
                univers = rest.split(" ")[-1]
    # build diamond database
    cmd = ["diamond", "makedb", "--in", "uniprot_sprot.fasta", "--db", "uniprot"]
    logger.info(f"Building diamond database: {' '.join(cmd)}")
    subprocess.call(cmd, cwd=os.path.join(env["FUNANNOTATE2_DB"]))
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
    # delete existintg
    for f in os.listdir(env["FUNANNOTATE2_DB"]):
        if f.startswith("dbCAN."):
            os.remove(os.path.join(env["FUNANNOTATE2_DB"], f))
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("dbCAN"))
    # download data
    download(DBURL.get("dbCAN"), hmm, wget=wget)
    download(DBURL.get("dbCAN-tsv"), familyinfo, wget=wget)
    download(DBURL.get("dbCAN-log"), versionfile, wget=wget)
    dbdate = ""
    dbvers = ""
    with io.open(versionfile, encoding="utf8", errors="ignore") as infile:
        head = [next(infile) for x in range(2)]
    dbdate = head[1].replace("# ", "").rstrip()
    dbvers = head[0].split(" ")[-1].rstrip()
    dbdate = datetime.datetime.strptime(dbdate, "%m/%d/%Y").strftime("%Y-%m-%d")
    logger.info("Creating dbCAN HMM database and pressing with pyhmmer")
    num_records = pyhmmer.hmmer.hmmpress(pyhmmer.plan7.HMMFile(hmm), hmm)
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
    for f in os.listdir(env["FUNANNOTATE2_DB"]):
        if f.startswith("Pfam"):
            os.remove(os.path.join(env["FUNANNOTATE2_DB"], f))
    logger.info("Downloading Pfam database")
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("pfam-log"))
    download(DBURL.get("pfam"), hmm + ".gz", wget=wget)
    download(DBURL.get("pfam-tsv"), familyinfo + ".gz", wget=wget)
    download(DBURL.get("pfam-log"), versionfile + ".gz", wget=wget)
    subprocess.call(["gunzip", "-f", "Pfam-A.hmm.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"]))
    subprocess.call(
        ["gunzip", "-f", "Pfam-A.clans.tsv.gz"],
        cwd=os.path.join(env["FUNANNOTATE2_DB"]),
    )
    subprocess.call(["gunzip", "-f", "Pfam.version.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"]))
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
    actual_recs = pyhmmer.hmmer.hmmpress(pyhmmer.plan7.HMMFile(hmm), hmm)
    assert actual_recs == num_records
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
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("go"))
    download(DBURL.get("go"), goOBO, wget=wget)
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
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("mibig"))
    download(DBURL.get("mibig"), fasta, wget=wget)
    version = os.path.basename(DBURL.get("mibig")).split("_")[-1].replace(".fasta", "")
    cmd = ["diamond", "makedb", "--in", "mibig.fa", "--db", "mibig"]
    logger.info(f"Building diamond database: {' '.join(cmd)}")
    subprocess.call(cmd, cwd=os.path.join(env["FUNANNOTATE2_DB"]))
    # runSubprocess(cmd, os.path.join(env['FUNANNOTATE2_DB']), logger)
    num_records = countfasta(fasta)
    return {
        "type": "diamond",
        "db": database,
        "version": version,
        "date": today,
        "n_records": num_records,
        "md5": md5,
    }


def mitoDB(wget=False):
    mitoFA = os.path.join(env["FUNANNOTATE2_DB"], "mito.refseq.fasta")
    mitoMMI = os.path.join(env["FUNANNOTATE2_DB"], "mito.mmi")
    mitoREF = os.path.join(env["FUNANNOTATE2_DB"], "mito.refseq.version")
    logger.info("Downloading NCBI RefSeq mitochondrial genomes")
    for x in [mitoFA, mitoMMI, mitoREF]:
        if checkfile(x):
            os.remove(x)
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("mito-release"))
    download(DBURL.get("mito"), mitoFA, wget=wget)
    download(DBURL.get("mito-release"), mitoREF, wget=wget)
    with io.open(mitoREF, encoding="utf8", errors="ignore") as infile:
        version = infile.readline().strip()
    # now create mmi database
    cmd = ["minimap2", "-x", "asm20", "-d", mitoMMI, mitoFA]
    subprocess.call(cmd, cwd=os.path.join(env["FUNANNOTATE2_DB"]))
    num_records = countfasta(mitoFA)
    # we can delete the fasta as not needed
    if checkfile(mitoFA):
        os.remove(mitoFA)
    return {
        "type": "minimap2",
        "db": mitoMMI,
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
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("interpro-tsv"))
    download(DBURL.get("interpro"), iprXML + ".gz", wget=wget)
    download(DBURL.get("interpro-tsv"), iprTSV, wget=wget)
    subprocess.call(["gunzip", "-f", "interpro.xml.gz"], cwd=os.path.join(env["FUNANNOTATE2_DB"]))
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
    # get md5 remote to store so update works properly
    md5 = calcmd5remote(DBURL.get("gene2product"))
    download(DBURL.get("gene2product"), curatedFile, wget=wget)
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
