# tools for storing and fetching training data
import json
import os
import shutil

from natsort import natsorted

from .config import env
from .log import startLogging, system_info
from .utilities import checkfile, print_table


def species(args):
    """
    Manage species data operations including loading, deleting, and displaying species information.

    This function handles various operations related to species data management in the FUNANNOTATE2_DB. It can load new species parameters into the database, delete existing species, and display the pre-trained species data. The function starts a logger, logs system information, and interacts with the database based on the provided command-line arguments.

    Parameters:
        args (argparse.Namespace): Command-line arguments parsed by argparse, containing options for loading, deleting, and formatting species data.

    Returns:
        None
    """
    # start logger
    logger = startLogging()
    system_info(logger.info)

    if args.load:  # then load new data
        logger.info(f"Loading new params file into database: {args.load}")
        loaded = load_params(args.load)
        if loaded is None:
            logger.error(
                "Species is not unique and already exists, to overwrite/delete try funannotate2 species --delete <name>"
            )
            raise SystemExit(1)
    if args.delete:  # then we are deleting a params/species
        existing = fetch_pretrained_species()
        if args.delete not in existing:
            logger.error(
                f"Species {args.delete} is not found in pre-trained database {os.path.join(env['FUNANNOTATE2_DB'], 'pretrained')}"
            )
            raise SystemExit(1)
        # now delete it
        del_params = os.path.join(
            os.path.join(env["FUNANNOTATE2_DB"], "pretrained", f"{args.delete}.params.json")
        )
        del_dir = os.path.join(os.path.join(env["FUNANNOTATE2_DB"], "pretrained", args.delete))
        if checkfile(del_params):
            os.remove(del_params)
        if os.path.isdir(del_dir):
            shutil.rmtree(del_dir)

    # at end always show the data in database
    db_species = show_species(layout=args.format)
    logger.info(f"Pre-trained species in database [format={args.format}]:\n{db_species}")
    logger.info(
        'Note, use these spcies in funannotate2 predict with the "-p, --pretrained" argument'
    )


def show_species(layout="table"):
    """
    Display species information in a specified layout.

    This function retrieves pretrained species data and displays it in either a table format or as a JSON string, depending on the specified layout.

    Parameters:
        layout (str, optional): The layout in which to display the species information. Defaults to "table". If "table", the information is formatted as a table; otherwise, it is returned as a JSON string.

    Returns:
        str: A formatted table of species information if `layout` is "table"; otherwise, a JSON string of all species.
    """
    all_species = fetch_pretrained_species()
    if layout == "table":
        d = [["name", "species", "busco-lineage", "abinitio-models", "date", "user"]]
        for k, v in natsorted(all_species.items()):
            bl = v["busco-lineage"]
            d.append(
                [
                    v["name"],
                    v["species"],
                    os.path.basename(bl),
                    ",".join(list(v["abinitio"].keys())),
                    v["date"].split()[0],
                    v["user"],
                ]
            )
        return print_table(d, return_str=True, max_col_width=40)
    else:
        return json.dumps(all_species, indent=2)


def fetch_pretrained_species():
    """
    Load existing or installed species parameters from the pretrained folder.

    This function accesses the `FUNANNOTATE2_DB` directory to retrieve species parameters stored in JSON files within the "pretrained" folder. It compiles these parameters into a dictionary, where each key is the species name.

    Returns:
        dict: A dictionary containing the loaded species parameters, with species names as keys.
    """
    # load existing/installed species/parameters
    existing = {}
    folder = os.path.join(env["FUNANNOTATE2_DB"], "pretrained")
    if os.path.isdir(folder):
        for f in os.listdir(folder):
            if f.endswith(".params.json"):
                with open(os.path.join(folder, f)) as infile:
                    load = json.load(infile)
                    if load["name"] not in existing:
                        existing[load["name"]] = load
    return existing


def load_params(paramfile, existing={}):
    """
    Load parameters from a params.json file and install them in the FUNANNOTATE2_DB directory.

    This function reads a JSON file containing species parameters and checks if the species is already present in the existing parameters. If the species is new, it copies necessary files to the appropriate directory and updates the parameters.

    Parameters:
        paramfile (str): The path to the params.json file containing species parameters.
        existing (dict, optional): A dictionary containing loaded species parameters. Defaults to an empty dictionary.

    Returns:
        str: The path where the updated parameters are saved if the species is new, otherwise None.
    """
    # takes a params.json file and installs in $FUNANNOTATE2_DB
    if not existing:
        existing = fetch_pretrained_species()
    with open(paramfile, "r") as infile:
        params = json.load(infile)
    if params["name"] not in existing:  # then its new
        updated = {
            "name": params["name"],
            "species": params["species"],
            "taxonomy": params["taxonomy"],
            "busco-lineage": params["busco-lineage"],
            "date": params["date"],
            "user": params["user"],
            "abinitio": {},
        }
        # directory to save data to and update pointers
        sp_dir = os.path.join(env["FUNANNOTATE2_DB"], "pretrained", params["name"])
        if not os.path.isdir(sp_dir):
            os.makedirs(sp_dir)
        for k, v in params["abinitio"].items():
            if k == "augustus":
                new_loc = os.path.join(sp_dir, "augustus")
                shutil.copytree(v["location"], os.path.abspath(new_loc))
                n = v.copy()
                n["location"] = os.path.abspath(new_loc)
                n.pop("training_set", None)
                updated["abinitio"]["augustus"] = n
            elif k == "glimmerhmm":
                new_loc = os.path.join(sp_dir, "glimmerhmm")
                shutil.copytree(v["location"], os.path.abspath(new_loc))
                n = v.copy()
                n["location"] = os.path.abspath(new_loc)
                n.pop("training_set", None)
                updated["abinitio"]["glimmerhmm"] = n
            elif k == "snap":
                new_loc = os.path.join(sp_dir, "snap", "snap-trained.hmm")
                if not os.path.isdir(os.path.join(sp_dir, "snap")):
                    os.makedirs(os.path.join(sp_dir, "snap"))
                shutil.copyfile(v["location"], os.path.abspath(new_loc))
                n = v.copy()
                n["location"] = os.path.abspath(new_loc)
                n.pop("training_set", None)
                updated["abinitio"]["snap"] = n
            elif k == "genemark":
                new_loc = os.path.join(sp_dir, "genemark", "gmhmm.mod")
                if not os.path.isdir(os.path.join(sp_dir, "genemark")):
                    os.makedirs(os.path.join(sp_dir, "genemark"))
                shutil.copyfile(v["location"], os.path.abspath(new_loc))
                n = v.copy()
                n["location"] = os.path.abspath(new_loc)
                n.pop("training_set", None)
                updated["abinitio"]["genemark"] = n
        # now can write the json file
        final_place = os.path.join(
            env["FUNANNOTATE2_DB"], "pretrained", os.path.basename(paramfile)
        )
        with open(final_place, "w") as outfile:
            json.dump(updated, outfile, indent=2)
        return final_place
    else:
        return None
